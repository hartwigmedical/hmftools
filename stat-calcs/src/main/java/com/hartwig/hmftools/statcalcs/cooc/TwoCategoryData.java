package com.hartwig.hmftools.statcalcs.cooc;

public class TwoCategoryData
{
    public final String Category1;
    public final String Category2;
    public int Count;

    public TwoCategoryData(final String category1, final String category2, int count)
    {
        Category1 = category1;
        Category2 = category2;
        Count = count;
    }
}
