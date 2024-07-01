package com.hartwig.hmftools.pave;

import java.util.List;

import com.google.common.collect.Lists;

public class PaveConstants
{
    public static final String APP_NAME = "Pave";

    public static final int GENE_UPSTREAM_DISTANCE = 1000;
    public static final int SPLICE_REGION_INTRON_RANGE = 8;
    public static final int SPLICE_REGION_EXON_RANGE = 3;

    public static final int SPLICE_ACCEPTOR_END_RANGE = 3;
    public static final List<Integer> SPLICE_ACCEPTOR_POSITIONS = Lists.newArrayList(1, 2, SPLICE_ACCEPTOR_END_RANGE);
    public static final List<Integer> SPLICE_DONOR_POSITIONS = Lists.newArrayList(-1, 1, 2, 5);

    // currently only TERT sets coding impact upstream of the coding start
    public static final List<String> PROMOTOR_UPSTREAM_GENE_IDS = Lists.newArrayList("ENSG00000164362");
    public static final int PROMOTOR_UPSTREAM_DISTANCE = 300;

    // PON related
    public static final int PON_SAMPLE_COUNT_THRESHOLD = 10;
    public static final int PON_REPEAT_COUNT_THRESHOLD = 4;
    public static final int PON_MEAN_READ_THRESHOLD = 6;
    public static final double PON_VAF_THRESHOLD = 0.08;

    public static final double GNMOAD_FILTER_HOTSPOT_PATHOGENIC_THRESHOLD = 0.01;
    public static final double GNMOAD_FILTER_THRESHOLD = 0.00015;

    public static final int PON_FILTER_HOTSPOT_SAMPLE_COUNT = 10;
    public static final int PON_FILTER_HOTSPOT_MAX_READS = 5;

    public static final int PON_FILTER_PANEL_SAMPLE_COUNT = 6;
    public static final int PON_FILTER_PANEL_MAX_READS = 5;

    public static final int PON_FILTER_OTHER_TIER_SAMPLE_COUNT = 6;
    public static final int PON_FILTER_OTHER_TIER_MAX_READS = 0;
}
