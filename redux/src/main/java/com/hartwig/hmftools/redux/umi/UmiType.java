package com.hartwig.hmftools.redux.umi;

public enum UmiType
{
    NONE,
    SINGLE,
    TWIST_DUPLEX,
    TSO500_DUPEX,
    OTHER_DUPLEX;

    public static final String TSO500_DUPEX_DELIM = "+";
    public static final String TWIST_DUPEX_DELIM = "_";
}
