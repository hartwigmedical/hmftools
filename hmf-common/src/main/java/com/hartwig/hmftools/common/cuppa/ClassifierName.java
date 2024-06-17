package com.hartwig.hmftools.common.cuppa;

public enum ClassifierName
{
    COMBINED,

    DNA_COMBINED,
    GEN_POS,
    SNV96,
    EVENT,

    RNA_COMBINED,
    GENE_EXP,
    ALT_SJ,

    NONE;

    public static String convertAliasToName(String string)
    {
        string = string.toUpperCase();

        if(string.equals("RMD"))
        {
            return "GEN_POS";
        }

        return string;
    }

}
