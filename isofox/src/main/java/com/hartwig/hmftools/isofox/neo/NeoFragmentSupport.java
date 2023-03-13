package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;

public class NeoFragmentSupport
{
    public static final int EXACT_MATCH = 0;
    public static final int PARTIAL_MATCH = 1;
    public static final int MISMATCH = 2;

    public int[] UpFragments;
    public int[] NovelFragments;
    public int[] DownFragments;

    public int[] RefBaseDepth;

    public NeoFragmentSupport()
    {
        RefBaseDepth = new int[] { 0, 0 };
        UpFragments = new int[] { 0, 0 };
        NovelFragments = new int[] { 0, 0 };
        DownFragments = new int[] { 0, 0 };
    }

    public boolean hasAnySupport()
    {
        return hasSupport(EXACT_MATCH) || hasSupport(PARTIAL_MATCH);
    }

    public boolean hasSupport(int type)
    {
        return UpFragments[type] > 0 || NovelFragments[type] > 0 || DownFragments[type] > 0;
    }

    public void combineSupport(final NeoFragmentSupport other)
    {
        for(int i = EXACT_MATCH; i <= PARTIAL_MATCH; ++i)
        {
            UpFragments[i] += other.UpFragments[i];
            NovelFragments[i] += other.NovelFragments[i];
            DownFragments[i] += other.DownFragments[i];
        }
    }

    public void setMax(final NeoFragmentSupport other)
    {
        for(int i = EXACT_MATCH; i <= PARTIAL_MATCH; ++i)
        {
            UpFragments[i] = max(UpFragments[i], other.UpFragments[i]);
            NovelFragments[i] = max(NovelFragments[i], other.NovelFragments[i]);
            DownFragments[i] = max(DownFragments[i], other.DownFragments[i]);
        }
    }

    public String toString()
    {
        return format("depth(%d:%d) up(%d:%d) novel(%d:%d) novel(%d:%d)",
                RefBaseDepth[FS_UP], RefBaseDepth[FS_DOWN], UpFragments[FS_UP], UpFragments[FS_DOWN],
                NovelFragments[FS_UP], NovelFragments[FS_DOWN], DownFragments[FS_UP], DownFragments[FS_DOWN]);
    }
}
