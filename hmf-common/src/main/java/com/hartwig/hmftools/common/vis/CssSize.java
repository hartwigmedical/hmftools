package com.hartwig.hmftools.common.vis;

import static java.lang.String.format;

public class CssSize
{
    public static final CssSize ZERO = px(0);

    public final double Magnitude;
    public final CssUnit Unit;

    private CssSize(double magnitude, final CssUnit unit)
    {
        Magnitude = magnitude;
        Unit = unit;
    }

    public static CssSize px(double magnitude)
    {
        return new CssSize(magnitude, CssUnit.PX);
    }

    public static CssSize em(double magnitude)
    {
        return new CssSize(magnitude, CssUnit.EM);
    }

    public boolean isZero()
    {
        return Math.round(10000 * Magnitude) == 0;
    }

    public CssSize scale(double factor)
    {
        return new CssSize(factor * Magnitude, Unit);
    }

    @Override
    public String toString()
    {
        if(isZero())
        {
            return "0";
        }

        return format("%.4f", Magnitude) + Unit;
    }

    private enum CssUnit
    {
        PX("px"),
        EM("em");

        private final String mString;

        private CssUnit(final String string)
        {
            mString = string;
        }

        @Override
        public String toString()
        {
            return mString;
        }
    }
}
