package com.hartwig.hmftools.common.codon;

public final class HgvsCommon
{
    public static final String HGVS_TYPE_DEL = "del";
    public static final String HGVS_TYPE_DUP = "dup";
    public static final String HGVS_TYPE_INS = "ins";

    public static final String HGVS_TYPE_DEL_INS = "delins";

    // coding
    public static final String HGVS_CODING_ID = "c.";
    public static final String HGVS_NON_CODING_ID = "n.";

    // protein
    public static final String HGVS_PROTEIN_ID = "p.";
    public static final String HGVS_FRAMESHIFT = "fs";
    public static final String HGVS_STOP_LOST = "ext*?";
    public static final String HGVS_START_LOST = "?";
    public static final String HGVS_STOP_GAINED = "*";
    public static final String HGVS_STOP_TRI_CODE = "Ter";
    public static final String HGVS_SYNONYMOUS = "=";
}
