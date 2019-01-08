package com.hartwig.hmftools.svannotation.analysis;

import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

public class RnaFusionData
{
    public final String Name;
    public final String GeneUp;
    public final String GeneDown;
    public final String ChrUp;
    public final String ChrDown;
    public final long PositionUp;
    public final long PositionDown;
    public final byte StrandUp;
    public final byte StrandDown;

    public final int JunctionReadCount;
    public final int SpanningFragCount;
    public final String SpliceType;

    private int mExonMinRankUp;
    private int mExonMaxRankUp;
    private int mExonMinRankDown;
    private int mExonMaxRankDown;

    private GeneFusion mFusionMatch;
    private Transcript mTransUp;
    private Transcript mTransDown;


    public RnaFusionData(final String name, final String geneUp, final String geneDown, final String chrUp, final String chrDown,
            long posUp, long posDown, byte strandUp, byte strandDown, int junctionReadCount, int spanningFragCount, final String spliceType)
    {
        Name = name;
        GeneUp = geneUp;
        GeneDown = geneDown;
        ChrUp = chrUp;
        ChrDown = chrDown;
        PositionUp = posUp;
        PositionDown = posDown;
        StrandUp = strandUp;
        StrandDown = strandDown;
        JunctionReadCount = junctionReadCount;
        SpanningFragCount = spanningFragCount;
        SpliceType = spliceType;

        mExonMinRankUp = 0;
        mExonMaxRankUp = 0;
        mExonMinRankDown = 0;
        mExonMaxRankDown = 0;

        mFusionMatch = null;
        mTransUp = null;
        mTransDown = null;
    }

    public void setExonUpRank(int min, int max)
    {
        mExonMaxRankUp = max;
        mExonMinRankUp = min;
    }

    public void setExonDownRank(int min, int max)
    {
        mExonMaxRankDown = max;
        mExonMinRankDown = min;
    }

    public int exonMinRankUp() { return mExonMinRankUp; }
    public int exonMaxRankUp() {return mExonMaxRankUp; }
    public int exonMinRankDown() { return mExonMinRankDown; }
    public int exonMaxRankDown() { return mExonMaxRankDown; }

    public final GeneFusion getMatchedFusion() { return mFusionMatch; }
    public void setMatchedFusion(final GeneFusion fusion)
    {
        mFusionMatch = fusion;
        mTransUp = fusion.upstreamTrans();
        mTransDown = fusion.downstreamTrans();
    }

    public void setBreakends(final Transcript up, final Transcript down)
    {
        mTransUp = up;
        mTransDown = down;
    }

    public final Transcript getTransUp() { return mTransUp; }
    public final Transcript getTransDown() { return mTransDown; }

}
