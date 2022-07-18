package com.hartwig.hmftools.pave;

import java.util.List;

import com.google.common.collect.Lists;

public class PaveConstants
{
    public static final int GENE_UPSTREAM_DISTANCE = 1000;
    public static final int SPLICE_REGION_INTRON_RANGE = 8;
    public static final int SPLICE_REGION_EXON_RANGE = 3;

    public static final int SPLICE_ACCEPTOR_END_RANGE = 3;
    public static final List<Integer> SPLICE_ACCEPTOR_POSITIONS = Lists.newArrayList(1, 2, SPLICE_ACCEPTOR_END_RANGE);
    public static final List<Integer> SPLICE_DONOR_POSITIONS = Lists.newArrayList(-1, 1, 2, 5);

    // currently only TERT sets coding impact upstream of the coding start
    public static final List<String> PROMOTOR_UPSTREAM_GENE_IDS = Lists.newArrayList("ENSG00000164362");
    public static final int PROMOTOR_UPSTREAM_DISTANCE = 300;

    public static final String DELIM = ",";
    public static final String ITEM_DELIM = ";";
}
