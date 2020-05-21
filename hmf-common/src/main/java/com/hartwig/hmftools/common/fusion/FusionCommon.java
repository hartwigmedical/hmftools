package com.hartwig.hmftools.common.fusion;

public class FusionCommon
{
    public static final int FS_UPSTREAM = 0;
    public static final int FS_DOWNSTREAM = 1;
    public static final int FS_PAIR = 2;

    public static final String UPSTREAM_STR = "upstream";
    public static final String DOWNSTREAM_STR = "downstream";

    public static int switchStream(int iter) { return iter == FS_UPSTREAM ? FS_DOWNSTREAM : FS_UPSTREAM; }
    public static String streamStr(int iter) { return iter == FS_UPSTREAM ? UPSTREAM_STR : DOWNSTREAM_STR; }

}
