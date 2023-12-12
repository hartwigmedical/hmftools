import numpy as np
import pandas as pd
from sklearn.feature_selection import chi2
import scipy as sp
from scipy.stats import ranksums, median_abs_deviation

def chi2_test(X, y, sort_by_stat=True):
    results = chi2(X, y)

    output = pd.DataFrame({
        "feature": X.columns,
        "stat": results[0],
        "pvalue": results[1],
        "rank": (-results[0]).argsort().argsort()
    })

    if sort_by_stat:
        output = output.sort_values("stat", ascending=False)

    return output


def wilcox_test_ovr(X, y, sort_by_stat=True, top_n=None):
    if False:
        top_n = 5
        sort_by_stat = True

    unique_classes = np.unique(y)
    output = []
    for class_i in unique_classes:
        # class_i=unique_classes[0]
        X_true = X.loc[y == class_i]
        X_false = X.loc[y != class_i]
        results = ranksums(X_true, X_false)
        output.append(
            pd.DataFrame({
                "class": class_i,
                "feature": X.columns,
                "stat": results[0],
                "pvalue": results[1],
                "rank": (-results[0]).argsort().argsort()
            })
        )

    output = pd.concat(output)

    if sort_by_stat:
        output = output.sort_values(["class", "stat"], ascending=[True, False])

    if top_n is not None:
        output = output[output["rank"]<=top_n]

    return output

def cohen_d(X1, X2, avg_metric="mean", sort_by_stat=True, top_n=None):
    if False:
        # M = get_X_columns(training_data, prefixes="gen_pos", return_data=True)
        # X1 = M.loc[(y=="Melanoma").values,:]
        # X2 = M.loc[(y!="Melanoma").values,:]

        X1 = X.loc[(y=="Lung: Non-small Cell").values,:]
        X2 = X.loc[(y!="Lung: Non-small Cell").values,:]

    ## https://stats.stackexchange.com/questions/497295/effect-size-calculation-for-comparison-between-medians
    ## (median(X1) - median(X2)) / sqrt((mad(X1)^2 + mad(X1)^2)/2) ## when X1 and X2 have same sample size

    if avg_metric=="mean":
        average_func = np.mean
        deviation_func = np.std
    elif avg_metric=="median":
        average_func = np.median
        deviation_func = median_abs_deviation
    else:
        raise Exception("`avg_metric` must be 'mean' or 'median'")

    X1_average = average_func(X1, axis=0)
    X2_average = average_func(X2, axis=0)

    X1_deviation = deviation_func(X1, axis=0)
    X2_deviation = deviation_func(X2, axis=0)

    ## Pooled (standard) deviation formula https://youtu.be/IetVSlrndpI?t=288
    n_samples_1 = X1.shape[0]
    n_samples_2 = X2.shape[0]

    stat = np.sqrt(
        (n_samples_1 - 1) * X1_deviation ** 2 + \
        (n_samples_2 - 1) * X2_deviation ** 2 / \

        (n_samples_1 + n_samples_2 - 2)
    )

    output = pd.DataFrame({
        "feature": X1.columns,
        "stat": stat,
        "rank": (-stat).argsort().argsort()

    }, index=None)

    if sort_by_stat:
        output = output.sort_values("stat", ascending=False)

    if top_n is not None:
        output = output[output["rank"] <= top_n]

    return output

def cohen_d_ovr(X, y, **kwargs):
    if False:
        top_n = 5
        sort_by_stat = True

    unique_classes = np.unique(y)
    output = []
    for class_i in unique_classes:
        # class_i=unique_classes[0]
        X_true = X.loc[y == class_i]
        X_false = X.loc[y != class_i]
        results = cohen_d(X_true, X_false)
        results.insert(0, "class", class_i)
        output.append(results)

    output = pd.concat(output)

    return output

# df = cohen_d_ovr(X, y)
#
# cohen_d(X.loc[(y=="Melanoma").values,:], X.loc[(y!="Melanoma").values,:], avg_metric="median", top_n=5)
# cohen_d(X.loc[(y=="Lung: Non-small Cell").values,:], X.loc[(y!="Lung: Non-small Cell").values,:])
# cohen_d(X.loc[(y=="Lymphoid tissue").values,:], X.loc[(y!="Lymphoid tissue").values,:])
# df = X.describe().T

def cliff_delta_array(x1, x2):
    if False:
        x1 = np.array(sigs.loc[y=="Melanoma","SBS7a"])
        x2 = np.array(sigs.loc[y!="Melanoma","SBS7a"])

    sign_sums = np.sum(np.sign(np.subtract.outer(x1, x2)))
    n_comparisons = x1.size * x2.size
    stat = sign_sums / n_comparisons

    return stat

def cliff_delta_matrix(X1, X2):
    if False:
        M = get_X_columns(training_data, prefixes="gen_pos", return_data=True)
        X1 = M.loc[(y=="Melanoma").values,:]
        X2 = M.loc[(y!="Melanoma").values,:]

    ## Rows: features. Columns: samples
    X1 = np.array(X1).T
    X2 = np.array(X2).T

    ##
    if X1.shape[0] != X2.shape[0]:
        raise Exception("Input matrices must have the same number of features (columns)")

    n_features = X1.shape[0]
    sign_sums = np.zeros(n_features)

    for i in range(0, n_features):
        #i=0
        sign_sum = np.sum(np.sign(np.subtract.outer(X1[i], X2[i])))
        sign_sums[i] = sign_sum

    n_comparisons = X1.shape[1] * X2.shape[1]
    stat = sign_sums / n_comparisons
    return stat


def fisher_test(case_true=None, case_false=None, ctrl_true=None, ctrl_false=None, matrix=None, alternative='two-sided'):
    """Perform a vectorized fishers exact test

    Parameters
    ------
    case_true, case_false, ctrl_true, ctrl_false: array-like
        Each argument is an array representing a position
        in the below 2x2 contingency matrix layout:
            case_true   ctrl_true
            case_false  ctrl_false
    matrix: numpy 2 array or pandas data frame
        Alternatively, the above arguments can be provided as
        a matrix with the form where each row is a sample:
            case_true_1  case_false_1  ctrl_true_1  ctrl_false_1
            case_true_2  case_false_2  ctrl_true_2  ctrl_false_2
            ...

    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis. Default is 'two-sided'.
        The following options are available:
        * 'two-sided': case is statistically greater or less than ctrl
        * 'less': case is statistically less than ctrl
        * 'greater': case is statistically greater than ctrl

    Returns
    -------
    An array of pvalues

    Examples
    --------
    ##
    case_true=[6,100,5]
    case_false=[1,10,5]
    ctrl_true=[2,2,5]
    ctrl_false=[4,69,5]
    fisher_test(
        case_true, case_false, ctrl_true, ctrl_false
    )
    (0.10256410256410256, 2.764512960449913e-36, 1.0)

    Description
    ------
    This function is a reimplementation of the code from the
    docs of `scipy.stats.fisher_exact`:
        table = np.array([[6, 2], [1, 4]]) ## np.array([[case_true, case_false], [ctrl_true, ctrl_false]])
        M = table.sum()
        n = table[0].sum()
        N = table[:, 0].sum()
        start, end = hypergeom.support(M, n, N)
        hypergeom.pmf(np.arange(start, end+1), M, n, N)

    """
    ## Example input
    # if False:
    #     case_true = [6, 100, 5]
    #     case_false = [1, 10, 5]
    #     ctrl_true = [2, 2, 5]
    #     ctrl_false = [4, 69, 5]
    #     matrix = np.column_stack(
    #         (case_true, case_false, ctrl_true, ctrl_false)
    #     )
    #     matrix = pd.DataFrame(matrix)
    #     alternative = 'two-sided'

    ## Convert matrix to 4 arrays ----------------------------
    if matrix is not None:
        if matrix.shape[1] != 4:
            raise ValueError(
                '`matrix` must have 4 columns corresponding to `case_true`, `case_false`, `ctrl_true`, and `ctrl_false`')

        case_true, case_false, ctrl_true, ctrl_false = matrix[0], matrix[1], matrix[2], matrix[3]

    ## Check inputs ----------------------------
    ## Convert inputs to int array
    case_true, case_false, ctrl_true, ctrl_false = map(
        np.asarray,
        (case_true, case_false, ctrl_true, ctrl_false)
    )

    ## Check if inputs have the same length
    vector_lengths = (case_true, case_false, ctrl_true, ctrl_false)
    if len(set(i.size for i in vector_lengths)) != 1:
        raise ValueError('The 4 input vectors must be of the same length')
    del vector_lengths

    ## Alternative hypothesis
    n_tests = len(case_true)
    if isinstance(alternative, str):
        alternative = [alternative] * n_tests
    elif len(alternative) != n_tests:
        raise ValueError("`alternative` must be the same length as the input data")

    if not all(i in ('two-sided', 'greater', 'less') for i in alternative):
        raise ValueError("`alternative` must be 'two-sided', 'greater', or 'less'")

    ## Main ----------------------------
    M = case_true + case_false + ctrl_true + ctrl_false
    n = case_true + case_false
    N = case_true + ctrl_true

    def calc_fisher_test_pvalue(M, n, N, case_true, alternative):
        start, end = sp.stats.hypergeom.support(M, n, N)
        x = np.arange(start, end + 1)
        p = sp.stats.hypergeom.pmf(x, M, n, N)

        if alternative == 'two-sided':
            p_thres = p[x == case_true]
            return sum(p[p <= p_thres])

        if alternative == 'greater':
            return sum(p[x >= case_true])

        if alternative == 'less':
            return sum(p[x <= case_true])

    # calc_fisher_test_pvalue(M[0], n[0], N[0], case_true[0], alternative[0])

    return tuple(map(calc_fisher_test_pvalue, M, n, N, case_true, alternative))

def contingency_table_ovr(X, y):
    # X_fusion = X[make_column_selector("fusion")].copy()
    # X = X_fusion

    unique_classes = np.unique(y)
    output = []
    for class_i in unique_classes:
        X_true = X.loc[y == class_i]
        X_false = X.loc[y != class_i]

        case_true = X_true.sum(axis=0)
        case_false = X_true.shape[0] - case_true

        ctrl_true = X_false.sum(axis=0)
        ctrl_false = X_false.shape[0] - case_false

        tab = pd.DataFrame({
            "class": class_i,
            "feature": X.columns,
            "case_true": case_true,
            "case_false": case_false,
            "ctrl_true": ctrl_true,
            "ctrl_false": ctrl_false
        })
        tab = tab.reset_index(drop=True)
        output.append(tab)

    return pd.concat(output)

def fisher_test_ovr(X, y, alternative="two-sided"):
    conting = contingency_table_ovr(X,y)
    conting["pvalue"] = fisher_test(
        case_true=conting["case_true"],
        case_false=conting["case_false"],
        ctrl_true=conting["ctrl_true"],
        ctrl_false=conting["ctrl_false"],
        alternative=alternative
    )
    return conting