#!/usr/bin/env python3

import argparse
import re
import sys
from abc import ABC, abstractmethod

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import stats

SCHEMA = pa.schema([("variant", pa.string()), ("beta", pa.float64()),
                    ("standard_error", pa.float64()), ("i_squared", pa.float64()),
                    ("sample_size", pa.float64())])


class QtlFilter(ABC):
    _operator = None
    _field_name = None

    @abstractmethod
    def _get_value(self):
        raise NotImplementedError

    def get_filter(self):
        return self._field_name, self._operator, self._get_value()

    def get_field_name(self):
        return self._field_name


class QtlSetFilter(QtlFilter):
    _operator = "in"

    def __init__(self, values):
        self._values = set(values)

    def _get_value(self):
        return self._values

    def get_values(self):
        return self._get_value()

    @classmethod
    def from_path(cls, path):
        df = pd.read_csv(path, header = None, delimiter = '\t')
        return cls(df.iloc[:,0].tolist())

    @classmethod
    def from_list(cls, genes):
        return cls(genes)

    def apply_to_df(self, df):
        return df.loc[df[self._field_name].isin(self._values)]


class QtlVariantFilter(QtlSetFilter):
    _field_name = "variant"
    def __init__(self, variants):
        super().__init__(variants)


class QtlGeneFilter(QtlSetFilter):
    _field_name = "phenotype"
    def __init__(self, genes):
        super().__init__(genes)


class QtlChromosomeFilter(QtlFilter):
    _field_name = "chromosome"
    _operator = "=="
    def __init__(self, chromosome):
        self._value = chromosome

    def _get_value(self):
        return self._value


class QtlLocusVariantFilter:

    def __init__(self, chromosome, variants):
        self._chromosome_filter = QtlChromosomeFilter(chromosome)
        self._variant_filter = QtlVariantFilter(variants)

    def get_filters(self):
        # Return a list, with in each list a chromosome filter, and a
        return [self._chromosome_filter, self._variant_filter]


class QtlPThresholdFilter(QtlFilter):
    _field_name = "p_value"
    _operator = "<"

    def __init__(self, p_value):
        self._value = p_value

    def apply_to_df(self, df):
        return df.loc[df[self._field_name] < self._value]

    def apply(self, series):
        return series < self._value

    def _get_value(self):
        return self._value


class QtlResultProcessor:

    def __init__(self, path, gene_filter=None):
        self.path = path
        self.gene_filter = gene_filter
        self.variant_filters = None
        self.significance_filter = None
        self.column_mapping = \
            {"t_stat": self.t_stat,
             "z_score": self.z_score,
             "p_value": self.p_value
             }
        self._t_stat = None
        self._z_score = None
        self._p_value = None
        self._df = None

    def extract(self,
                cols=None,
                drop=None,
                add=None):

        filters = self.get_filters()

        print("Using filters: {}".format(filters))

        dataset = pq.ParquetDataset(
            self.path, validate_schema=True,
            filters=filters)

        self._df = dataset.read().to_pandas()
        print("Queried {} rows".format(self._df.shape[0]))

        if self.significance_filter is not None:
            self._df = self._df.loc[self.significance_filter.apply(self.p_value())]
            self._p_value = self._p_value[self.significance_filter.apply(self.p_value())]

        print("Filtered to {} rows".format(self._df.shape[0]))

        cols_to_add = (
            cols.union(add)
            .intersection(self.column_mapping.keys()))

        print("Columns to add: {}".format(cols_to_add))

        self.include_cols(cols_to_add)

        print(self._df.head())

        print("Added columns. Table now has {} columns".format(self._df.shape[1]))

        if cols is not None and len(cols) > 0:
            out = self._df.loc[:, ["variant", "phenotype"].extend(list(cols))]
        elif drop is not None and len(drop) > 0:
            out = self._df.drop(drop, axis=1)
        else:
            out = self._df
        return out

    def t_stat(self):
        if self._t_stat is None:
            self._t_stat = self.calculate_t_stat()
        return self._t_stat

    def p_value(self):
        if self._p_value is None:
            self._p_value = self.calculate_p_value()
        return self._p_value

    def z_score(self):
        if self._z_score is None:
            self._z_score = self.calculate_t_stat()
        return self._z_score

    def calculate_t_stat(self):
        return self._df['beta']/self._df['standard_error']

    def calculate_z_score(self):
        return self.calculate_t_stat()

    def calculate_p_value(self, degrees_of_freedom = None):
        if degrees_of_freedom is None:
            return stats.norm.sf(np.abs(self.z_score()))*2
        else:
            return stats.t.sf(np.abs(self.z_score()), degrees_of_freedom)*2

    def get_filters(self):
        base_filter_list = [self.gene_filter.get_filter()]
        if self.variant_filters is not None:
            filter_list = [base_filter_list.extend(var_filter.get_filters()) for var_filter in self.variant_filters]
        else:
            filter_list = [base_filter_list]
        return filter_list

    def include_cols(self, cols):
        for col in cols:
            self._df[col] = self.column_mapping[col]()

def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description = "Analyze HASE output parquet files.")

    parser.add_argument('-i', '--input-file', required = True,
                        help = """One or multiple input parquet arrays. 
                                  Usage of wildcards is supported but then the argument has to be quoted.""")
    parser.add_argument('-o', '--output-file', type = str,
                        required = True,
                        help = "Path to tab-separated output file.")
    parser.add_argument('-p', '--p-thresh', type = float,
                        required = False,
                        help = "P-value threshold for filtering the results. Defaults None.",
                        default = None)
    parser.add_argument('-g', '--genes', required = False, default = None, nargs = '+',
                        help = """Individual phenotype IDs specified and separated by space.""")
    parser.add_argument('-G', '--genes-file', required = False, default = None, type = argparse.FileType('r'),
                        help = """File with the list of phenotypes to include.""")
    parser.add_argument('-v', '--variants', required = False, default = None, nargs = '+',
                        help = """Individual SNP IDs specified and separated by space.""")
    parser.add_argument('-V', '--variants-file', required = False, default = None, type = argparse.FileType('r'),
                        help = """File with the list of SNPs/variants to include.""")
    parser.add_argument('-r', '--variant-reference', required = False, type = argparse.FileType('r'),
                        help = "Reference for variants. Has to be gzipped and space-delimited.")
    parser.add_argument('-c', '--cols', required = False, default = None,
                        help = """Extract only z-scores""")

    args = parser.parse_args(argv[1:])

    print(args)

    variant_reference = None
    variants_list = None
    variant_filters = None

    if args.variant_reference is not None:
        variant_reference = (
            pd.read_csv(args.variant_reference, compression = 'gzip', sep = ' ')
            .drop(["allele1", "allele2"])
            .rename({"ID": "variant", "bp": "bp", "CHR": "chromosome", "str_allele1": "a1", "str_allele2": "a2"}))

    if args.variants_file is not None:
        print("Using variants file '%s' to filter on variants." % args.variants_file)
        variants_list = (
            pd.read_csv(args.variants_file, header = None, delimiter = '\t')
            .iloc[:,0].tolist())
    if args.variants is not None:
        print("Provided %d variants for filtering." % len(args.variants))
        if variants_list is not None:
            print("Variant filter already defined. Skipping...")
        variants_list = args.variants

    if variants_list is not None:
        if variant_reference is None:
            parser.error("Cannot subset on variants without variant reference")
        variant_selection = variant_reference.loc[variant_reference.variant.isin(variants_list), :]
        variant_dictionary = variant_selection.groupby('chromosome')['variant'].apply(list).to_dict()
        variant_filters = [QtlLocusVariantFilter(chromosome, variants) for chromosome, variants in variant_dictionary.items()]

    qtl_gene_filter = None

    if args.genes_file is not None:
        print("Using variants file '%s' to filter on variants." % args.genes_file)
        qtl_gene_filter = QtlGeneFilter.from_path(args.genes_file)
    if args.genes is not None:
        print("Provided %d genes for filtering." % len(args.genes))
        if qtl_gene_filter is not None:
            print("Variant filter already defined. Skipping...")
        qtl_gene_filter = QtlGeneFilter.from_list(args.genes)

    column_specifications = {"-": set(), "+": set(), "": set()}
    if args.cols is not None and args.cols != "":
        for column_specification in args.cols.split(","):
            regex_match = re.match(r"([+-]?)(\w+)", column_specification)
            if not regex_match:
                parser.error(
                    "Column specification {} did not match expected pattern"
                    .format(column_specification))
            column_specifications[regex_match.group(1)].add(regex_match.group(2))

    print(column_specifications)

    first = True

    with open(args.output_file, 'w') as f:

        for gene in qtl_gene_filter.get_values():
            print("Gene {}".format(gene))

            qtl_single_gene_filter = QtlGeneFilter.from_list([gene])

            result_processor = QtlResultProcessor(
                args.input_file, qtl_single_gene_filter)
            result_processor.variant_filters = variant_filters
            result_processor.significance_filter = QtlPThresholdFilter(args.p_thresh)

            df = result_processor.extract(
                cols=column_specifications[""],
                drop=column_specifications["-"],
                add=column_specifications["+"])
            print("out:")
            print(df.shape)
            print(df.head())
            df.to_csv(f, sep="\t", header=first, index=None)

            first = False

        print("Done!")
        print("Closing output file '{}'".format(args.output_file))

    return 0


if __name__ == '__main__':
    sys.exit(main())
