package com.hartwig.hmftools.linx.rna;

import static com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion.REPORTABLE_TYPE_NONE;
import static com.hartwig.hmftools.linx.types.SvChain.CHAIN_LENGTH;
import static com.hartwig.hmftools.linx.types.SvChain.CHAIN_LINK_COUNT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

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
    public static String RNA_SPLICE_TYPE_UNKONWN = "UNKNOWN"; // for read with a star-fusion prediction

    private boolean mIsValid;
    private List<Integer> mPotentialTransUp; // IDs of transcripts with an exon matching the RNA position
    private List<Integer> mPotentialTransDown;

    // annotations and matching results

    private Transcript mTransUp;
    private Transcript mTransDown;

    // canonical exon positions
    private int mExonMinRankUp;
    private int mExonMaxRankUp;
    private int mExonMinRankDown;
    private int mExonMaxRankDown;

    private boolean mViableFusion; // the pair of transcripts satisfied standard fusion rules
    private boolean mPhaseMatchedFusion;
    private String mKnownFusionType;
    private boolean mTransViableUp; // the transcript fell in the correct location relative to the RNA position
    private boolean mTransViableDown;
    private boolean mTransCorrectLocationUp; // the transcript fell on the correct side and orientation for the RNA position
    private boolean mTransCorrectLocationDown;
    private int mExonsSkippedUp; // where no valid breakend was found, record the number of exons skipped between the breakend and the RNA
    private int mExonsSkippedDown;

    // SVA match data
    private SvBreakend mBreakendUp;
    private SvBreakend mBreakendDown;

    private String mClusterInfoUp;
    private String mClusterInfoDown;
    private String mChainInfo;

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

        mIsValid = true;
        mExonMinRankUp = 0;
        mExonMaxRankUp = 0;
        mExonMinRankDown = 0;
        mExonMaxRankDown = 0;

        mPotentialTransDown = Lists.newArrayList();
        mPotentialTransUp = Lists.newArrayList();

        mTransUp = null;
        mTransDown = null;
        mBreakendUp = null;
        mBreakendDown = null;

        mViableFusion = false;
        mPhaseMatchedFusion = false;
        mKnownFusionType = REPORTABLE_TYPE_NONE;
        mTransViableUp = false;
        mTransViableDown = false;
        mTransCorrectLocationUp = false;
        mTransCorrectLocationDown = false;
        mExonsSkippedUp = 0;
        mExonsSkippedDown = 0;

        mClusterInfoUp = "";
        mClusterInfoDown = "";
        mChainInfo = "0;0";
    }

    public boolean isValid() { return mIsValid; }
    public void setValid(boolean toggle) { mIsValid = toggle; }

    public void setExonRank(boolean isUpstream, int min, int max)
    {
        if(isUpstream)
        {
            mExonMaxRankUp = max;
            mExonMinRankUp = min;
        }
        else
        {
            mExonMaxRankDown = max;
            mExonMinRankDown = min;
        }
    }

    public void setPotentialTrans(boolean isUpstream, final List<Integer> list)
    {
        if(isUpstream)
            mPotentialTransUp.addAll(list);
        else
            mPotentialTransDown.addAll(list);
    }

    public final List<Integer> getPotentialTrans(boolean isUpstream) { return isUpstream ? mPotentialTransUp : mPotentialTransDown; }

    public int exonMinRankUp() { return mExonMinRankUp; }
    public int exonMaxRankUp() {return mExonMaxRankUp; }
    public int exonMinRankDown() { return mExonMinRankDown; }
    public int exonMaxRankDown() { return mExonMaxRankDown; }

    public void setViableFusion(boolean viable, boolean phaseMatched)
    {
        mViableFusion = viable;
        mPhaseMatchedFusion = phaseMatched;
    }

    public void setTranscriptData(boolean isUpstream, final Transcript trans, final SvBreakend breakend,
            boolean matchedRnaBoundary, boolean correctLocation, int exonsSkipped)
    {
        if(isUpstream)
        {
            mTransUp = trans;
            mBreakendUp = breakend;
            mTransViableUp = matchedRnaBoundary;
            mTransCorrectLocationUp = correctLocation;
            mExonsSkippedUp = exonsSkipped;
        }
        else
        {
            mTransDown = trans;
            mBreakendDown = breakend;
            mTransViableDown = matchedRnaBoundary;
            mTransCorrectLocationDown = correctLocation;
            mExonsSkippedDown = exonsSkipped;
        }
    }

    public final Transcript getTrans(boolean isUpstream) { return isUpstream ? mTransUp : mTransDown; }
    public final SvBreakend getBreakend(boolean isUpstream) { return isUpstream ? mBreakendUp : mBreakendDown; }
    public final boolean isTransViable(boolean isUpstream) { return isUpstream ? mTransViableUp : mTransViableDown; }
    public final boolean isTransCorrectLocation(boolean isUpstream) { return isUpstream ? mTransCorrectLocationUp : mTransCorrectLocationDown; }
    public final int getExonsSkipped(boolean isUpstream) { return isUpstream ? mExonsSkippedUp : mExonsSkippedDown; }

    public boolean isViableFusion() { return mViableFusion; }
    public boolean isPhaseMatchedFusion() { return mPhaseMatchedFusion; }

    public void setKnownType(final String knownType) { mKnownFusionType = knownType; }
    public final String getKnownType() { return mKnownFusionType; }

    public final String getClusterInfo(boolean isUpstream)
    {
        return isUpstream ? mClusterInfoUp : mClusterInfoDown;
    }

    public final String getChainInfo()
    {
        return mChainInfo;
    }

    public void setFusionClusterChainInfo()
    {
        if(mBreakendUp == null && mBreakendDown == null)
            return;

        if(mBreakendUp != null && mBreakendDown != null)
        {
            SvVarData varUp = mBreakendUp.getSV();
            SvCluster clusterUp = varUp.getCluster();
            SvVarData varDown = mBreakendDown.getSV();
            SvCluster clusterDown = varDown.getCluster();

            SvChain matchingChain = null;

            if(varUp != varDown && clusterUp == clusterDown)
            {
                // check for a matching chain if the clusters are the same
                matchingChain = clusterUp.findSameChainForSVs(varUp, varDown);

                if (matchingChain != null)
                {
                    final int chainData[] =
                            matchingChain.breakendsAreChained(varUp, !mBreakendUp.usesStart(), varDown, !mBreakendDown.usesStart());
                    mChainInfo = String.format("%d;%d", chainData[CHAIN_LINK_COUNT], chainData[CHAIN_LENGTH]);
                }
            }

            setFusionClusterInfo(mBreakendUp, true, matchingChain);
            setFusionClusterInfo(mBreakendDown, false, matchingChain);
        }
        else
        {
            SvBreakend breakend = mBreakendUp != null ? mBreakendUp : mBreakendDown;
            setFusionClusterInfo(breakend, mBreakendUp != null, null);
        }
    }

    private void setFusionClusterInfo(final SvBreakend breakend, boolean isUpstream, SvChain matchingChain)
    {
        // data: ClusterId;ClusterCount;ChainId;ChainCount
        SvCluster cluster = breakend.getCluster();

        SvChain chain = matchingChain != null ? matchingChain : cluster.findChain(breakend.getSV());

        final String clusterInfo = String.format("%d;%d;%d;%d",
                cluster.id(), cluster.getSvCount(),
                chain != null ? chain.id() : cluster.getChainId(breakend.getSV()), chain != null ? chain.getSvCount() : 1);

        if(isUpstream)
            mClusterInfoUp = clusterInfo;
        else
            mClusterInfoDown = clusterInfo;
    }

}
