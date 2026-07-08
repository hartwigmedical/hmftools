package com.hartwig.hmftools.compar.common;

public class FieldOverride
{
    public final String Category;
    public final String Field;
    public final String Compared;
    public final String AbsoluteThreshold;
    public final String PercentThreshold;

    public FieldOverride(
            final String category, final String field, final String compared,
            final String absoluteThreshold, final String percentThreshold)
    {
        Category = category;
        Field = field;
        Compared = compared;
        AbsoluteThreshold = absoluteThreshold;
        PercentThreshold = percentThreshold;
    }
}
