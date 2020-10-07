package com.hartwig.hmftools.isofox.common;

public enum FragmentMatchType
{
    SPLICED,
    SHORT,
    LONG,
    UNSPLICED,
    DISCORDANT,
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
            case DISCORDANT: return 4;
            default: return 5;
        }
    }
}
