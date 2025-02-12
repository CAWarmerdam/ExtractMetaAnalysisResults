import numpy as np
import pandas as pd
from math import floor
from scipy.optimize import minimize_scalar
from itertools import combinations
from scipy.linalg import inv, det, pinv


def CARMA_spike_slab_noEM(z, ld, lambda_val=1, Max_Model_Dim=2e5,
                          all_iter=3, all_inner_iter=10, epsilon_threshold=1e-5, num_causal=10,
                          y_var=1, tau=0.04, outlier_switch=True, outlier_BF_index=1/3.2):

    """
    Perform CARMA analysis using a Spike-and-Slab prior without Expectation-Maximization (EM).

    Parameters:
    - z (array-like): Numeric vector representing z-scores.
    - ld (array-like): Numeric matrix representing the linkage disequilibrium (LD) matrix.
    - lambda_val (float): Regularization parameter controlling the strength of the L1 penalty.
    - Max_Model_Dim (float): Maximum allowed dimension for the causal models.
    - all_iter (int): The total number of iterations to run the CARMA analysis.
    - all_inner_iter (int): The number of inner iterations in each CARMA iteration.
    - epsilon_threshold (float): Threshold for convergence in CARMA iterations.
    - num_causal (int): Maximal number of causal variants to be selected in the final model.
    - y_var (float): Variance parameter for the marginal distribution of the response variable.
    - tau (float): Tuning parameter controlling the degree of sparsity in the Spike-and-Slab prior.
    - outlier_switch (bool): Whether to consider outlier detection in the analysis.
    - outlier_BF_index (float): Bayes Factor threshold for identifying outliers.

    Returns:
    - dict: A dictionary containing the following components:
        - 'PIPs': Vector of Posterior Inclusion Probabilities (PIPs) for each genetic variant.
        - 'B_list': Data frame containing information about identified models/configurations.
        - 'Outliers': List of identified outliers.
    """

    t0 = pd.Timestamp.now()

    p_snp = len(z)
    epsilon_list = epsilon_threshold * p_snp
    all_epsilon_threshold = epsilon_threshold * p_snp

    print("ordering pizza...")
    print(f'N SNPs is {p_snp}')
    # Zero step
    all_C_list = MCS_modified(z=z, ld_matrix=ld, epsilon=epsilon_list,
                              Max_Model_Dim=Max_Model_Dim, lambda_val=lambda_val,
                              outlier_switch=outlier_switch, tau=tau,
                              num_causal=num_causal, y_var=y_var,
                              inner_all_iter=all_inner_iter, outlier_BF_index=outlier_BF_index)

    t1 = pd.Timestamp.now() - t0
    print(f'Zero step is finished in {round(t1.total_seconds(), 2)} sec')
    print(f'Expecting to finish in {round(t1.total_seconds()*all_iter, 2)} sec')

    outliers=all_C_list["conditional_S_list"]
    # Main steps
    for g in range(0, all_iter):

        ac1 = all_C_list["B_list"]["set_gamma_margin"]
        previous_result = np.mean(ac1[0:round(len(ac1) / 4)])

        all_C_list = MCS_modified(z=z, ld_matrix=ld, input_conditional_S_list=all_C_list["conditional_S_list"],
                                  Max_Model_Dim=Max_Model_Dim,
                                  y_var=y_var, num_causal=num_causal, epsilon=epsilon_list,
                                  outlier_switch=outlier_switch, tau=tau,
                                  lambda_val=lambda_val,
                                  inner_all_iter=all_inner_iter,
                                  outlier_BF_index=outlier_BF_index)

        ac1 = all_C_list["B_list"]["set_gamma_margin"]
        difference = np.abs(previous_result - np.mean(ac1[0:round(len(ac1) / 4)]))
        if difference < all_epsilon_threshold:
            break

    # Calculate PIPs and Credible Set
    pip = PIP_func(likeli=all_C_list["B_list"]["set_gamma_margin"], model_space=all_C_list["B_list"]["matrix_gamma"], p=p_snp,num_causal=num_causal)

    results_list = {"PIPs":pip, "B_list":all_C_list["B_list"], "Outliers":all_C_list["conditional_S_list"]}

    print("pizza ordered!")
    t1 = pd.Timestamp.now() - t0

    print(f'FM time is {round(t1.total_seconds(), 2)} sec')

    return results_list



def ind_normal_sigma_fixed_marginal_fun_indi(zSigmaz_S, tau, p_S, det_S):
    result = p_S / 2.0 * np.log(tau) - 0.5 * np.log(det_S) + zSigmaz_S / 2.0
    return result

def ind_Normal_fixed_sigma_marginal_external(index_vec_input, Sigma, z, tau, p_S, y_sigma):
    index_vec = index_vec_input - 1
    Sigma_S = Sigma[np.ix_(index_vec, index_vec)]
    A = tau * np.eye(p_S)

    det_S = det(Sigma_S + A)
    Sigma_S_inv = inv(Sigma_S + A)

    sub_z = z[index_vec]
    zSigmaz_S = np.dot(sub_z.T, np.dot(Sigma_S_inv, sub_z))

    b = ind_normal_sigma_fixed_marginal_fun_indi(zSigmaz_S, tau, p_S, det_S)

    results = b

    return results


def outlier_ind_Normal_marginal_external(index_vec_input, Sigma, z, tau, p_S, y_sigma):
    index_vec = index_vec_input - 1

    Sigma_S = Sigma[np.ix_(index_vec, index_vec)]
    A = tau * np.eye(p_S)

    Sigma_S_I_inv = pinv(Sigma_S + A, rtol=0.00001)
    Sigma_S_inv = pinv(Sigma_S, rtol=0.00001)

    det_S = np.abs(det(Sigma_S_inv))
    det_I_S = np.abs(det(Sigma_S_I_inv))

    sub_z = z[index_vec]
    zSigmaz_S = np.dot(sub_z, np.dot(Sigma_S_inv, sub_z))
    zSigmaz_I_S = np.dot(sub_z, np.dot(Sigma_S_I_inv, sub_z))

    b = 0.5 * (np.log(det_S) + np.log(det_I_S)) - 0.5 * (zSigmaz_S - zSigmaz_I_S)
    results = b

    return results

def add_function(S_sub, y):
    results = [np.sort(np.concatenate(([x],y))) for x in S_sub]
    return np.array(results)

def set_gamma_func_base(S, p):
    set_gamma = {}

    # set of gamma-
    if len(S) == 0:
        set_gamma[0] = None
        set_gamma[1] = np.arange(0, p ).reshape(-1, 1)
        set_gamma[2] = None

    if len(S) == 1:
        S_sub = np.setdiff1d(np.arange(0, p ), S)
        set_gamma[0] = None
        set_gamma[1] = add_function(S_sub, S)
        set_gamma[2] = S_sub.reshape(-1, 1)

    if len(S) > 1:
        S_sub = np.setdiff1d(np.arange(0, p ), S)
        S = np.sort(S)
        set_gamma[0] = np.array(list(combinations(S, len(S) - 1)))
        set_gamma[1] = add_function(S_sub, S)
        xs = np.vstack([add_function(S_sub, row) for row in set_gamma[0]])
        set_gamma[2] = xs

    return set_gamma

def set_gamma_func_conditional(input_S, condition_index, p):
    set_gamma = {}
    S = np.setdiff1d(input_S, condition_index)

    # set of gamma-
    if len(S) == 0:
        S_sub = np.setdiff1d(np.arange(0, p), condition_index)
        set_gamma[0] = None
        set_gamma[1] = S_sub.reshape(-1, 1)
        set_gamma[2] = None

    if len(S) == 1:
        S_sub = np.setdiff1d(np.arange(0, p), input_S)
        set_gamma[0] = None
        set_gamma[1] = add_function(S_sub, S)
        set_gamma[2] = S_sub.reshape(-1, 1)

    if len(S) > 1:
        S_sub = np.setdiff1d(np.arange(0, p), input_S)
        S = np.sort(S)
        set_gamma[0] = np.array(list(combinations(S, len(S) - 1)))
        set_gamma[1] = add_function(S_sub, S)
        xs = np.vstack([add_function(S_sub, row) for row in set_gamma[0]])
        set_gamma[2] = xs

    return set_gamma

def set_gamma_func(input_S, p, condition_index=None):
    if condition_index is None:
        results = set_gamma_func_base(input_S, p)
    else:
        results = set_gamma_func_conditional(input_S, condition_index, p)
    return results

def index_fun_internal(x):
    y=np.sort(x)
    y=y.astype(str)
    l=','.join(y)
    return l

def index_fun(y):
    l=np.array([index_fun_internal(x) for x in y])
    return l

def ridge_fun(x, Sigma, modi_ld_S, test_S, z, outlier_tau,outlier_likelihood):
    temp_Sigma = Sigma.copy()
    temp_ld_S = x * modi_ld_S + (1 - x) * np.eye(len(modi_ld_S))
    temp_Sigma[np.ix_(test_S, test_S)] = temp_ld_S
    l = outlier_likelihood(index_vec_input=test_S + 1, Sigma=temp_Sigma, z=z, tau=outlier_tau, p_S=len(test_S), y_sigma=1)
    return -l

def prior_dist(t, lambda_val, p):
    l = t.split(",")
    dim_model = len(l)
    if t=="": dim_model=0
    result = dim_model * np.log(lambda_val) + np.math.lgamma(p - dim_model+1) - np.math.lgamma(p+1)
    return result

def PIP_func(likeli, model_space, p, num_causal):

    likeli=likeli.reset_index(drop=True)
    model_space=model_space.reset_index(drop=True)

    model_space_matrix = np.zeros((len(model_space), p), dtype=int)

    for i in range(len(model_space)):
        if model_space.iloc[i]!="":
            ind = list(map(int, model_space.iloc[i].split(',')))
            if len(ind) > 0:
                model_space_matrix[i, ind] = 1

    infi_index = np.where(np.isinf(likeli))[0]
    if len(infi_index) != 0:
        likeli=likeli.drop(infi_index).reset_index(drop=True)
        model_space_matrix = np.delete(model_space_matrix, infi_index, axis=0)

    na_index = np.where(np.isnan(likeli))[0]
    if len(na_index) != 0:
        likeli = likeli.drop(na_index).reset_index(drop=True)
        model_space_matrix = np.delete(model_space_matrix, na_index, axis=0)

    row_sums = np.sum(model_space_matrix, axis=1)
    model_space_matrix=model_space_matrix[row_sums<=num_causal]
    likeli=likeli[row_sums<=num_causal]

    aa = likeli - max(likeli)
    prob_sum = np.sum(np.exp(aa))


    result_prob = np.zeros(p)
    for i in range(p):
        result_prob[i] = np.sum(np.exp(aa[model_space_matrix[:, i] == 1])) / prob_sum

    return result_prob


def MCS_modified(z, ld_matrix, Max_Model_Dim=1e+4, lambda_val=1,
                 num_causal=10, y_var=1, outlier_switch=True,
                 input_conditional_S_list=None, tau=1/0.05**2,
                 epsilon=1e-3, inner_all_iter=10, outlier_BF_index=None):

    p=len(z)
    marginal_likelihood = ind_Normal_fixed_sigma_marginal_external
    tau_sample = tau
    if (outlier_switch == True):
        outlier_likelihood = outlier_ind_Normal_marginal_external
        outlier_tau = tau

    B = Max_Model_Dim
    stored_bf = 0
    Sigma =ld_matrix

    S = []

    null_model = ""
    null_margin = prior_dist(null_model, lambda_val=lambda_val, p=p)

    B_list = pd.DataFrame(
        {'set_gamma_margin': [null_margin],
         'matrix_gamma': [""]})

    if input_conditional_S_list==None:
        conditional_S=[]
    else:
        conditional_S = input_conditional_S_list
        S = conditional_S

    h_ind=0
    ind_inner_iter=0
    for ind_inner_iter in range(0,inner_all_iter):
        for h_ind in range(0,10):
            set_gamma=set_gamma_func(input_S=S, p=p, condition_index=conditional_S)

            if conditional_S is None:
                working_S=S
            else:
                working_S = np.sort(np.setdiff1d(S, conditional_S)).astype(int)

            set_gamma_margin = {}
            set_gamma_prior = {}
            matrix_gamma = {}

            if len(working_S)!=0:
                S_model = ','.join(np.sort(working_S).astype(int).astype(str))
                p_S = len(working_S)
                current_log_margin = marginal_likelihood(np.array(working_S)+1, Sigma, z, tau=tau_sample, p_S=p_S, y_sigma=y_var)+prior_dist(S_model, lambda_val=lambda_val, p=p)
            else:
                current_log_margin = prior_dist(null_model, lambda_val=lambda_val, p=p)

            set_gamma_margin = [None,None,None]
            set_gamma_prior = [None,None,None]
            matrix_gamma = [None,None,None]

            for i in range(0, len(set_gamma)):
                if (set_gamma[i] is not None):
                    matrix_gamma[i]=index_fun(set_gamma[i])
                    p_S = set_gamma[i].shape[1]
                    set_gamma_margin[i]=np.apply_along_axis(marginal_likelihood, 1, set_gamma[i]+1, Sigma=Sigma, z=z, tau=tau_sample, p_S=p_S,
                                                            y_sigma=y_var)
                    set_gamma_prior[i]=np.array([prior_dist(model, lambda_val=lambda_val, p=p) for model in matrix_gamma[i]])
                    set_gamma_margin[i] = set_gamma_prior[i] + set_gamma_margin[i]
                else:
                    set_gamma_margin[i]=np.array(null_margin)
                    set_gamma_prior[i]=0
                    matrix_gamma[i]=np.array(null_model)

            columns = ['set_gamma_margin', 'matrix_gamma']
            add_B = pd.DataFrame(columns=columns)

            for i in range(len(set_gamma)):
                if type(set_gamma_margin[i].tolist())==list:
                    new_row = pd.DataFrame({'set_gamma_margin': set_gamma_margin[i].tolist(),
                                            'matrix_gamma': matrix_gamma[i].tolist()})
                    add_B = pd.concat([add_B, new_row], ignore_index=True)
                else:
                    new_row = pd.DataFrame({'set_gamma_margin': [set_gamma_margin[i].tolist()],
                                            'matrix_gamma': [matrix_gamma[i].tolist()]})
                    add_B = pd.concat([add_B, new_row], ignore_index=True)


            # Add visited models into the storage space of models
            B_list = pd.concat([B_list, add_B], ignore_index=True)
            B_list=B_list.drop_duplicates(subset='matrix_gamma',ignore_index=True)
            B_list=B_list.sort_values(by='set_gamma_margin',ignore_index=True, ascending=False)

            if len(working_S)==0:
                # Create a DataFrame set.star
                set_star = pd.DataFrame({
                    'set_index': [0, 1, 2],
                    'gamma_set_index': [np.nan, np.nan, np.nan],
                    'margin': [np.nan, np.nan, np.nan]
                })

                # Assuming set.gamma.margin and current.log.margin are defined
                aa = set_gamma_margin[1]
                aa = aa - aa[np.argmax(aa)]

                min_half_len = min(len(aa), floor(p/2))
                decr_ind=np.argsort(np.exp(aa))[::-1]
                decr_half_ind=decr_ind[:min_half_len]

                probs=np.exp(aa)[decr_half_ind]

                chosen_index=np.random.choice(decr_half_ind, 1, p=probs/np.sum(probs))
                set_star.at[1, 'gamma_set_index'] = chosen_index[0]
                set_star.at[1, 'margin'] = set_gamma_margin[1][chosen_index[0]]

                S = set_gamma[1][chosen_index[0]].tolist()

            else:
                set_star = pd.DataFrame({
                    'set_index': [0, 1, 2],
                    'gamma_set_index': [np.nan, np.nan, np.nan],
                    'margin': [np.nan, np.nan, np.nan]
                })
                for i in range(0,3):
                    aa = set_gamma_margin[i]
                    if np.size(aa)>1:
                        aa = aa - aa[np.argmax(aa)]
                        chosen_index = np.random.choice(range(0, np.size(set_gamma_margin[i])), 1, p=np.exp(aa)/np.sum(np.exp(aa)))
                        set_star.at[i, 'gamma_set_index'] = chosen_index
                        set_star.at[i, 'margin'] = set_gamma_margin[i][chosen_index]
                    else:
                        chosen_index=0
                        set_star.at[i, 'gamma_set_index'] = chosen_index
                        set_star.at[i, 'margin'] = set_gamma_margin[i]

                if outlier_switch:
                    for i in range(1, len(set_gamma)):
                        test_log_BF=100
                        while True:
                            aa = set_gamma_margin[i]
                            aa = aa - aa[np.argmax(aa)]
                            chosen_index = np.random.choice(range(0, np.size(set_gamma_margin[i])), 1, p=np.exp(aa)/np.sum(np.exp(aa)))
                            set_star.at[i, 'gamma_set_index'] = chosen_index
                            set_star.at[i, 'margin'] = set_gamma_margin[i][chosen_index]

                            test_S = set_gamma[i][int(chosen_index),:]

                            modi_Sigma = Sigma.copy()
                            temp_Sigma = Sigma.copy()
                            if np.size(test_S) > 1:
                                modi_ld_S = modi_Sigma[test_S][:,test_S]

                                result =  minimize_scalar(ridge_fun, bounds=(0,1), args=(Sigma, modi_ld_S, test_S, z, outlier_tau,outlier_likelihood), method='bounded')
                                modi_ld_S = result.x * modi_ld_S + (1 - result.x) * np.eye(len(modi_ld_S))

                                modi_Sigma[np.ix_(test_S,test_S)] = modi_ld_S

                                test_log_BF=outlier_likelihood(test_S+1, Sigma, z, outlier_tau, len(test_S), 1) - outlier_likelihood(test_S+1, modi_Sigma, z, outlier_tau, len(test_S), 1)
                                test_log_BF = -np.abs(test_log_BF)

                            if np.exp(test_log_BF) < outlier_BF_index:
                                set_gamma[i] = np.delete(set_gamma[i], int(set_star['gamma_set_index'][i]), axis=0)
                                set_gamma_margin[i] = np.delete(set_gamma_margin[i], int(set_star['gamma_set_index'][i]), axis=0)
                                conditional_S = np.concatenate([conditional_S, np.setdiff1d(test_S, working_S)])
                                conditional_S = np.unique(conditional_S).astype(int).tolist()
                            else:
                                break


                if len(working_S) == num_causal:
                    set_star = set_star.drop(1)
                    aa = set_star['margin']- max(set_star['margin'])
                    sec_sample = np.random.choice([0, 2], 1, p=np.exp(aa) / np.sum(np.exp(aa)))
                    ind_sec=int(set_star["gamma_set_index"][set_star['set_index']==int(sec_sample)])
                    S=set_gamma[sec_sample[0]][ind_sec].tolist()
                else:
                    aa = set_star['margin']- max(set_star['margin'])
                    sec_sample = np.random.choice(range(0, 3), 1, p=np.exp(aa) / np.sum(np.exp(aa)))
                    S=set_gamma[sec_sample[0]][int(set_star['gamma_set_index'][sec_sample[0]])].tolist()

            for item in conditional_S:
                if item not in S:
                    S.append(item)
        #END h_ind loop
        #
        if conditional_S is not None:
            all_c_index = []
            l = [s.split(",") for s in B_list["matrix_gamma"]]
            for tt in conditional_S:
                tt=str(tt)
                ind = [i for i, sublist in enumerate(l) if tt in sublist]
                all_c_index.extend(ind)

            all_c_index = list(set(all_c_index))

            if len(all_c_index) > 0:
                temp_B_list = B_list.copy()
                temp_B_list = B_list.drop(all_c_index)
            else:
                temp_B_list = B_list.copy()
        else:
            temp_B_list = B_list.copy()

        result_B_list=temp_B_list[:min(int(B), len(temp_B_list))]

        rb1 = result_B_list["set_gamma_margin"]

        difference = abs(rb1[:(len(rb1)//4)].mean()-stored_bf)

        if difference < epsilon:
            break
        else:
            stored_bf = rb1[:(len(rb1)//4)].mean()

    out = {'B_list': result_B_list, 'conditional_S_list': conditional_S}

    return out