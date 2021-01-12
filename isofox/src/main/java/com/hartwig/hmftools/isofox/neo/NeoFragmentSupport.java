package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;

public class NeoFragmentSupport
{
    public static final int EXACT_MATCH = 0;
    public static final int PARTIAL_MATCH = 1;
    public static final int MISMATCH = 2;

    public int[] UpFragments;
    public int[] NovelFragments;
    public int[] DownFragments;

    public int[] RefBaseSupport;

    public NeoFragmentSupport()
    {
        RefBaseSupport = new int[] { 0, 0 };
        UpFragments = new int[] { 0, 0 };
        NovelFragments = new int[] { 0, 0 };
        DownFragments = new int[] { 0, 0 };
    }

    public void combine(final NeoFragmentSupport other)
    {
        for(int i = EXACT_MATCH; i <= PARTIAL_MATCH; ++i)
        {
            RefBaseSupport[i] += other.RefBaseSupport[i];
            UpFragments[i] += other.UpFragments[i];
            NovelFragments[i] += other.NovelFragments[i];
            DownFragments[i] += other.DownFragments[i];
        }
    }

    public void setMax(final NeoFragmentSupport other)
    {
        for(int i = EXACT_MATCH; i <= PARTIAL_MATCH; ++i)
        {
            RefBaseSupport[i] = max(RefBaseSupport[i], other.RefBaseSupport[i]);
            UpFragments[i] = max(UpFragments[i], other.UpFragments[i]);
            NovelFragments[i] = max(NovelFragments[i], other.NovelFragments[i]);
            DownFragments[i] = max(DownFragments[i], other.DownFragments[i]);
        }
    }

}
