package com.hartwig.hmftools.cup.somatics;

public final class SomaticsCommon
{
    public static final String SPLIT_AID_APOBEC = "split_aid_apobec_gen_pos";
    public static final String SPLIT_AID_APOBEC_DESC = "Exclude 8 AID/APOBEC trinucleotide contexts from genomic positions";

    public static final String NORMALISE_COPY_NUMBER = "normalise_cn";
    public static final String NORMALISE_COPY_NUMBER_DESC = "Adjust genomic-position counts by copy number ";

    public static final String EXCLUDE_SNV_96_AID_APOBEC = "exclude_aid_apobec_snv_96";
    public static final String EXCLUDE_SNV_96_AID_APOBEC_DESC = "Exclude 8 AID/APOBEC trinucleotide contexts from SNV-96 counts";

    public static final String INCLUDE_AID_APOBEC_SIG = "aid_apobec_sig_feature";
    public static final String INCLUDE_AID_APOBEC_SIG_DESC = "Add an enriched AID/APOBEC signature feature";
}
