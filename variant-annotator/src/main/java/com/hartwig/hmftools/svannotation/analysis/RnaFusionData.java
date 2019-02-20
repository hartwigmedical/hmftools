package com.hartwig.hmftools.svannotation.analysis;

import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

public class RnaFusionData
{
    // data from Star Fusion predictions output:
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

    public static String RNA_SPLICE_TYPE_ONLY_REF = "ONLY_REF_SPLICE";

    // annotations and matching results

    private Transcript mTransUp;
    private Transcript mTransDown;

    // canonical exon positions
    private int mExonMinRankUp;
    private int mExonMaxRankUp;
    private int mExonMinRankDown;
    private int mExonMaxRankDown;

    private boolean mViableFusion; // the pair of transcripts satisfied standard fusion rules
    private boolean mTransValidUp; // the transcript fell in the correct location relative to the RNA position
    private boolean mTransValidDown;

    // SVA match data
    private int mClusterIdUp;
    private int mClusterIdDown;
    private int mChainIdUp;
    private int mChainIdDown;

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

        mTransUp = null;
        mTransDown = null;

        mViableFusion = false;
        mTransValidUp = false;
        mTransValidDown = false;
        mClusterIdUp = -1;
        mClusterIdDown = -1;
        mChainIdUp = -1;
        mChainIdDown = -1;
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

    public void setViableFusion(boolean toggle) { mViableFusion = toggle; }

    public void setTranscriptData(final Transcript trans, boolean isUpstream, boolean matchedRnaBoundary)
    {
        if(isUpstream)
        {
            mTransUp = trans;
            mTransValidUp = matchedRnaBoundary;
        }
        else
        {
            mTransDown = trans;
            mTransValidDown = matchedRnaBoundary;
        }
    }

    public final Transcript getTrans(boolean isUpstream) { return isUpstream ? mTransUp : mTransDown; }
    public final boolean getTransValid(boolean isUpstream) { return isUpstream ? mTransValidUp : mTransValidDown; }

    public void setChainInfo(int clusterId, int chainId, boolean isUpstream)
    {
        if(isUpstream)
        {
            mChainIdUp = chainId;
            mClusterIdUp = clusterId;
        }
        else
        {
            mChainIdDown = chainId;
            mClusterIdDown = clusterId;
        }
    }

    public boolean isViableFusion() { return mViableFusion; }

    // SVA match data

    public int getClusterId(boolean isUpstream) { return isUpstream ? mClusterIdUp : mClusterIdDown; }
    public int getChainId(boolean isUpstream) { return isUpstream ? mChainIdUp : mChainIdDown; }

}
