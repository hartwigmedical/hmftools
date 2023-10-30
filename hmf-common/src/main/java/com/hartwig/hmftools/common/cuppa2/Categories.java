package com.hartwig.hmftools.common.cuppa2;

public class Categories
{
    public enum ClfName
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

        public static String convertAliasToName(String string){
            string = string.toUpperCase();

            if(string.equals("RMD")){
                return "GEN_POS";
            }

            return string;
        }

    }

    public enum ClfGroup
    {
        COMBINED,
        DNA,
        RNA,

        NONE;
    }

    public enum DataType
    {
        PROB,
        SIG_QUANTILE,
        FEAT_CONTRIB,
        CV_PERFORMANCE,

        NONE;
    }
}
