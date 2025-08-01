package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.redux.common.Constants.INVALID_BASE_QUAL;

public class BaseQualPair
{
    public final byte Base;
    public final byte Qual;

    public static final byte NO_BASE = 0;

    public static final BaseQualPair INVALID = new BaseQualPair(NO_BASE, INVALID_BASE_QUAL);

    public BaseQualPair(final byte base, final byte qual)
    {
        Base = base;
        Qual = qual;
    }

    public BaseQualPair(final byte base, final double qual)
    {
        this(base, (byte)round(qual));
    }

    public boolean isValid() { return Base == NO_BASE || Qual > 0; }

    public String toString() { return format("%d:%d", (char)Base, Qual); }
}
