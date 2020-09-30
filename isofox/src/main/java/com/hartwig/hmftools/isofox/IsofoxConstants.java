package com.hartwig.hmftools.isofox;

public class IsofoxConstants
{
    public static final int DEFAULT_MAX_FRAGMENT_SIZE = 550;

    // min number of fragments to sample when calculating fragment length distribution for calculation of expected fragment counts
    public static final int DEFAULT_FRAG_LENGTH_MIN_COUNT = 1000000;

    public static final double DEFAULT_GC_RATIO_BUCKET = 0.01;

    public static final short SINGLE_MAP_QUALITY = 255; // multi-mapped fragments are given map quals of 3 or lower

    public static final int ENRICHED_GENE_BUFFER = 100000;

    public static final int MAX_NOVEL_SJ_DISTANCE = 500000; // beyond which a fragment will be considered chimeric
}
