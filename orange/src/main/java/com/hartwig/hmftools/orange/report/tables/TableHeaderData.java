package com.hartwig.hmftools.orange.report.tables;

import static java.lang.Math.max;
import static java.lang.String.format;

import com.google.common.annotations.VisibleForTesting;

public class TableHeaderData
{
    public final String Text;
    public final double WrapAtPercent; // wrap the header text if required
    public final boolean CanWrapData; // can wrap the row data for this column if required

    public double HeaderWidth;
    public double DataWidth;
    public double ExcessWidth;

    public TableHeaderData(final String text, final double wrapAtPercent, final boolean canWrapData)
    {
        Text = text;
        WrapAtPercent = wrapAtPercent;
        CanWrapData = canWrapData;
        HeaderWidth = 0;
        DataWidth = 0;
        ExcessWidth = 0;
    }

    public int length() { return Text.length(); }

    public double adjustedWidth() { return max(HeaderWidth, DataWidth) + ExcessWidth; }

    public String toString()
    {
        return format("%s: wrapPerc(%.1f) dataWrap(%s) adjustedWidth(%.1f header=%.1f data=%.1f exess=%.1f)",
                Text, WrapAtPercent, CanWrapData, adjustedWidth(), HeaderWidth, DataWidth, ExcessWidth);
    }

    @VisibleForTesting
    public void reset()
    {
        HeaderWidth = 0;
        DataWidth = 0;
        ExcessWidth = 0;
    }
}
