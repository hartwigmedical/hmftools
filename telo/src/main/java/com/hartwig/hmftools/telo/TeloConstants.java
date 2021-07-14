package com.hartwig.hmftools.telo;

import java.util.List;

import com.google.common.collect.Lists;

public class TeloConstants
{
    public static final int DEFAULT_MIN_MAPPING_QUALITY = 10;
    public static final int DEFAULT_PARTITION_SIZE = 10000000;

    public static final int DEFAULT_MIN_TELE_SEQ_COUNT = 3;

    public static final String CANONICAL_TELOMERE_SEQ = "GGGTTA";
    public static final String CANONICAL_TELOMERE_SEQ_REV = "CCCTAA";

    public static final List<String> CANONICAL_TELOMERE_SEQUENCES = Lists.newArrayList(
            CANONICAL_TELOMERE_SEQ, CANONICAL_TELOMERE_SEQ_REV);

}
