package com.hartwig.hmftools.common.utils.sv;

public final class SvCommonUtils
{
    public static final byte POS_ORIENT = 1;
    public static final byte NEG_ORIENT = -1;

    public static final String POS_ORIENT_ID = "+";
    public static final String NEG_ORIENT_ID = "-";

    public static byte flipOrientation(final byte orientation) { return orientation == POS_ORIENT ? NEG_ORIENT : POS_ORIENT; }
}
