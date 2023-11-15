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

        public static String convertAliasToName(String string)
        {
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
        FEAT_CONTRIB,
        SIG_QUANTILE,
        CV_PERFORMANCE, // This data type is only use for plotting

        NONE;

        public static boolean isSampleLevelDataType(DataType dataType)
        {
            if(dataType.equals(PROB) | dataType.equals(FEAT_CONTRIB) | dataType.equals(SIG_QUANTILE))
            {
                return true;
            }
            return false;
        }
    }
}
