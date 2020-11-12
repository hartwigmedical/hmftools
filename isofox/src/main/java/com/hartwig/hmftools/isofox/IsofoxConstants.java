package com.hartwig.hmftools.isofox;

import com.hartwig.hmftools.common.utils.sv.SvRegion;

public class IsofoxConstants
{
    public static final int DEFAULT_MAX_FRAGMENT_SIZE = 550;

    // min number of fragments to sample when calculating fragment length distribution for calculation of expected fragment counts
    public static final int DEFAULT_FRAG_LENGTH_MIN_COUNT = 1000000;

    public static final double DEFAULT_GC_RATIO_BUCKET = 0.01;

    public static final short DEFAULT_SINGLE_MAP_QUALITY = 255;
    public static short SINGLE_MAP_QUALITY = DEFAULT_SINGLE_MAP_QUALITY;
    public static short MULTI_MAP_QUALITY_THRESHOLD = 3; // multi-mapped fragments are given map quals of 3 or lower

    public static final int ENRICHED_GENE_BUFFER = 100000;

    public static final int MAX_NOVEL_SJ_DISTANCE = 500000; // beyond which a fragment will be considered chimeric

    // LINC00486
    public static final SvRegion EXCLUDED_REGION_1_HG19 = new SvRegion("2", 33141260, 33141700);
    public static final SvRegion EXCLUDED_REGION_1_HG38 = new SvRegion("2", 32916190, 32916630);

    public static final String ENRICHED_GENE_1 = "ENSG00000258486"; // RN7SL1
    public static final String ENRICHED_GENE_2 = "ENSG00000265150"; // RN7SL2
    public static final String ENRICHED_GENE_3 = "ENSG00000266037"; // RN7SL3
    public static final String ENRICHED_GENE_4 = "ENSG00000263740"; // RN7SL4P
    public static final String ENRICHED_GENE_5 = "ENSG00000265735"; // RN7SL5P
    public static final String ENRICHED_GENE_6 = "ENSG00000202198"; // RN7SK
}
