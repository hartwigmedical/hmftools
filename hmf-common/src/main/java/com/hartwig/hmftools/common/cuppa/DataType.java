package com.hartwig.hmftools.common.cuppa;

public enum DataType
{
    PROB,
    FEAT_CONTRIB,
    SIG_QUANTILE,
    CV_PERFORMANCE, // This data type is only use for plotting

    NONE;

    public static boolean isSampleLevelDataType(DataType dataType)
    {
        return dataType.equals(PROB) || dataType.equals(FEAT_CONTRIB) || dataType.equals(SIG_QUANTILE);
    }
}
