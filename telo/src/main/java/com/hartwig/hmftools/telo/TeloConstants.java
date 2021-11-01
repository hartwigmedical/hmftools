package com.hartwig.hmftools.telo;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

public class TeloConstants
{
    public static final int DEFAULT_PARTITION_SIZE = 10000000;

    public static final int DEFAULT_MIN_TELE_SEQ_COUNT = 2;

    public static final double POLY_G_THRESHOLD = 0.9;

    public static final String CANONICAL_TELOMERE_SEQ = "TTAGGG";
    public static final String CANONICAL_TELOMERE_SEQ_REV = TeloUtils.reverseComplementSequence(CANONICAL_TELOMERE_SEQ);
    public static final String[] TELOMERE_HEXAMERS = {
            CANONICAL_TELOMERE_SEQ,
            "TCAGGG",
            "TTCGGG",
            "GTAGGG",
            "TGAGGG",
            "TTGGGG",
            "TAAGGG",
            "ATAGGG",
            "CTAGGG",
            "TTTGGG"
    };
    public static final String[] TELOMERE_HEXAMERS_REV = Arrays.stream(TELOMERE_HEXAMERS).map(TeloUtils::reverseComplementSequence).toArray(String[]::new);

    public static ChrBaseRegion[] EXCLUDED_BASE_REGIONS = {
            // LINC00486 ref 37
            new ChrBaseRegion("2", 33141260, 33141700),
            // LINC00486 ref 38
            new ChrBaseRegion("chr2", 32916190, 32916630)
    };
}
