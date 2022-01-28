package com.hartwig.hmftools.common.fusion;

public enum KnownFusionType
{
    NONE,
    KNOWN_PAIR,
    PROMISCUOUS_5,
    PROMISCUOUS_3,
    IG_KNOWN_PAIR,
    IG_PROMISCUOUS,
    EXON_DEL_DUP;

    public static final String PROMISCUOUS_BOTH = "PROMISCUOUS_BOTH";

    public static String displayStr(final String knownTypeStr)
    {
        if(knownTypeStr.equals(NONE.toString()))
            return "None";
        else if(knownTypeStr.equals(KNOWN_PAIR.toString()))
            return "Known pair";
        else if(knownTypeStr.equals(PROMISCUOUS_5.toString()))
            return "5' Promiscuous";
        else if(knownTypeStr.equals(PROMISCUOUS_3.toString()))
            return "3' Promiscuous";
        else if(knownTypeStr.equals(IG_KNOWN_PAIR.toString()))
            return "IG known pair";
        else if(knownTypeStr.equals(IG_PROMISCUOUS.toString()))
            return "IG promiscuous";
        else if(knownTypeStr.equals(EXON_DEL_DUP.toString()))
            return "Exon del dup";
        else if(knownTypeStr.equals(PROMISCUOUS_BOTH))
            return "5' and 3' Promiscuous";
        else
            return knownTypeStr;
    }
}
