package com.hartwig.hmftools.common.fusion;

public class FusionCommon
{
    public static final int DEFAULT_PRE_GENE_PROMOTOR_DISTANCE = 100000;

    public static final int FS_UP = 0;
    public static final int FS_DOWN = 1;
    public static final int FS_PAIR = 2;

    public static final byte POS_STRAND = 1;
    public static final byte NEG_STRAND = -1;

    public static final String UPSTREAM_STR = "upstream";
    public static final String DOWNSTREAM_STR = "downstream";

    public static int fsIndex(boolean isUpstream) { return isUpstream ? FS_UP : FS_DOWN; }

    public static int switchStream(int iter) { return iter == FS_UP ? FS_DOWN : FS_UP; }
    public static String streamStr(int iter) { return iter == FS_UP ? UPSTREAM_STR : DOWNSTREAM_STR; }

}
