package com.hartwig.hmftools.telo;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class TeloConstants
{
    public static final int DEFAULT_PARTITION_SIZE = 10000000;

    public static final int DEFAULT_MIN_TELE_SEQ_COUNT = 2;

    public static final String CANONICAL_TELOMERE_SEQ = "TTAGGG";
    public static final String CANONICAL_TELOMERE_SEQ_REV = "CCCTAA";
    public static final ChrBaseRegion UNMAPPED_BASE_REGION = new ChrBaseRegion("", 0, 0);
}
