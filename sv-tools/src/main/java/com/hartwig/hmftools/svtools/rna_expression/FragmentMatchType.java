package com.hartwig.hmftools.svtools.rna_expression;

public enum FragmentMatchType
{
    SPLICED,
    SHORT,
    LONG,
    UNSPLICED,
    MAX;

    public static final int MAX_FRAG_TYPE = typeAsInt(FragmentMatchType.MAX);

    public static int typeAsInt(FragmentMatchType type)
    {
        switch(type)
        {
            case SPLICED: return 0;
            case SHORT: return 1;
            case LONG: return 2;
            case UNSPLICED: return 3;
            default: return 4;
        }
    }

    public static FragmentMatchType intAsType(int type)
    {
        switch(type)
        {
            case 0: return SPLICED;
            case 1: return SHORT;
            case 2: return LONG;
            case 3: return UNSPLICED;
            default: return MAX;
        }
    }

}
