package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.ceil;
import static java.lang.String.format;

public class TableHeaderData
{
    public final String Text;
    public final double WrapAtPercent; // wrap the header text if required
    public final boolean CanWrapData; // can wrap the row data for this column if required

    public double AdjustedWidth;

    public TableHeaderData(final String text, final double wrapAtPercent, final boolean canWrapData)
    {
        Text = text;
        WrapAtPercent = wrapAtPercent;
        CanWrapData = canWrapData;
        AdjustedWidth = 0;
    }

    public int length() { return Text.length(); }
    public int minLength() { return WrapAtPercent == 0 ? length() : (int)ceil(WrapAtPercent * length()); }

    public String toString()
    {
        return format("%s: wrapPerc(%.1f) dataWrap(%s) adjustedWidth(%.1f)", Text, WrapAtPercent, CanWrapData, AdjustedWidth);
    }
}
