package com.hartwig.hmftools.linx.fusion.rna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_LENGTH;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.breakendsAreChained;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.KNOWN;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.utils.sv.StartEndPair;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class RnaFusionData
{
    // RNA fusion data from external source
    public final String SampleId;
    public final String FusionId;
    public final String Source;
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;
    public final byte[] Strands;

    public final int JunctionFragments;
    public final int DiscordantFragments;
    public final RnaJunctionType[] JunctionTypes;

    public final int CohortCount;
    public final String OtherData;

    // annotations
    private boolean mIsValid;

    // annotations and matching results

    // transcripts matching SV breakends
    private final BreakendTransData[] mMatchedTranscripts;

    private final StartEndPair<List<RnaExonMatchData>> mTransExonData;

    // transcripts matching RNA positions if in a phased fusion
    private final RnaExonMatchData[] mRnaPhasedTranscriptExons;

    private boolean mViableFusion; // the pair of transcripts satisfied standard fusion rules
    private boolean mPhaseMatchedFusion;
    private String mKnownFusionType;
    private final boolean[] mTransViable; // the transcript fell in the correct location relative to the RNA position
    private final boolean[] mTransCorrectLocation; // the transcript fell on the correct side and orientation for the RNA position
    private final int[] mExonsSkipped; // where no valid breakend was found, record the number of exons skipped between the breakend and the RNA

    // SVA match data
    private final SvBreakend[] mBreakends;

    private DnaRnaMatchType mDnaFusionMatchType; // no match, match on gene, match on exact SVs
    private String mDnaFusionMatchInfo;
    private boolean mDnaReportableFusion;

    private String mClusterInfoUp;
    private String mClusterInfoDown;
    private String mChainInfo;

    public RnaFusionData(
            final String sampleId, final String source, final String fusionId, final String[] geneIds, final String[] geneNames,
            final String[] chromosomes, final int[] positions, final byte[] orientations,
            int junctionReadCount, int spanningFragCount, final RnaJunctionType[] junctionTypes, int cohortCount, final String otherData)
    {
        SampleId = sampleId;
        Source = source;
        FusionId = fusionId;
        GeneIds = geneIds;
        GeneNames = geneNames;
        Chromosomes = chromosomes;
        Positions = positions;
        Orientations = orientations;
        JunctionFragments = junctionReadCount;
        DiscordantFragments = spanningFragCount;
        JunctionTypes = junctionTypes;
        CohortCount = cohortCount;
        OtherData = otherData;

        Strands = new byte[FS_PAIR];

        mIsValid = true;
        mRnaPhasedTranscriptExons = new RnaExonMatchData[] {null, null};

        mTransExonData = new StartEndPair<List<RnaExonMatchData>>( Lists.newArrayList(), Lists.newArrayList() );

        mMatchedTranscripts = new BreakendTransData[] { null, null};
        mBreakends = new SvBreakend[] { null, null};

        mViableFusion = false;
        mPhaseMatchedFusion = false;
        mKnownFusionType = KnownFusionType.NONE.toString();
        mTransViable = new boolean[] { false, false };
        mTransCorrectLocation = new boolean[] { false, false };
        mExonsSkipped = new int[] {0, 0};

        mDnaFusionMatchType = DnaRnaMatchType.NONE;
        mDnaFusionMatchInfo = "";
        mDnaReportableFusion = false;

        mClusterInfoUp = NO_CLUSTER_INFO;
        mClusterInfoDown = NO_CLUSTER_INFO;
        mChainInfo = NO_CHAIN_INFO;
    }

    public static final String NO_CLUSTER_INFO = "0;0;0;0";
    public static final String NO_CHAIN_INFO = "0;0";

    public String name() { return GeneNames[FS_UP] + "_" + GeneNames[FS_DOWN]; }

    public boolean isValid() { return mIsValid; }
    public void setValid(boolean toggle) { mIsValid = toggle; }

    public StructuralVariantType rnaSvType()
    {
        if(!Chromosomes[FS_UP].equals(Chromosomes[FS_DOWN]))
            return BND;
        if(Orientations[FS_UP] == Orientations[FS_DOWN])
            return INV;

        if(Positions[FS_UP] < Positions[FS_DOWN])
            return Orientations[FS_UP] == 1 ? DEL : DUP;
        else
            return Orientations[FS_UP] == -1 ? DEL : DUP;
    }

    public boolean matchesKnownSpliceSites() { return JunctionTypes[FS_UP] == KNOWN && JunctionTypes[FS_DOWN] == KNOWN; }

    public List<RnaExonMatchData> getTransExonData(int seIndex) { return mTransExonData.get(seIndex); }

    public final List<String> getExactMatchTransIds(int fs)
    {
        return mTransExonData.get(fs).stream().filter(x -> x.BoundaryMatch).map(x -> x.TransName).collect(Collectors.toList());
    }

    public RnaExonMatchData getBestExonMatch(int fs)
    {
        if(mRnaPhasedTranscriptExons[fs] != null)
            return mRnaPhasedTranscriptExons[fs];

        for(RnaExonMatchData exonData : mTransExonData.get(fs))
        {
            if(exonData.BoundaryMatch)
                return exonData;
        }

        return mTransExonData.get(fs).isEmpty() ? null : mTransExonData.get(fs).get(0);
    }

    public void setRnaPhasedFusionData(final RnaExonMatchData up, final RnaExonMatchData down)
    {
        mRnaPhasedTranscriptExons[FS_UP] = up;
        mRnaPhasedTranscriptExons[FS_DOWN] = down;
    }

    public boolean hasRnaPhasedFusion()
    {
        return mRnaPhasedTranscriptExons[FS_UP] != null && mRnaPhasedTranscriptExons[FS_DOWN] != null;
    }

    public void setViableFusion(boolean viable, boolean phaseMatched)
    {
        mViableFusion = viable;
        mPhaseMatchedFusion = phaseMatched;
    }

    public void setDnaFusionMatch(final DnaRnaMatchType matchType, final String matchInfo, boolean reportableFusion)
    {
        mDnaFusionMatchType = matchType;
        mDnaFusionMatchInfo = matchInfo;
        mDnaReportableFusion = reportableFusion;
    }

    public final DnaRnaMatchType getDnaFusionMatchType() { return mDnaFusionMatchType; }
    public final boolean hasDnaReportableFusion() { return mDnaReportableFusion; }

    public final String getDnaFusionMatchInfo()
    {
        return mDnaFusionMatchInfo.isEmpty() ?
                mDnaFusionMatchType.toString() : String.format("%s_%s", mDnaFusionMatchType, mDnaFusionMatchInfo);
    }

    public void setTranscriptData(int fs, final BreakendTransData trans, final SvBreakend breakend,
            boolean matchedRnaBoundary, boolean correctLocation, int exonsSkipped)
    {
        mMatchedTranscripts[fs] = trans;
        mBreakends[fs] = breakend;
        mTransViable[fs] = matchedRnaBoundary;
        mTransCorrectLocation[fs] = correctLocation;
        mExonsSkipped[fs] = exonsSkipped;
    }

    public final BreakendTransData[] getMatchedTranscripts() { return mMatchedTranscripts; }
    public final SvBreakend[] getBreakend() { return mBreakends; }
    public final boolean[] isTransViable() { return mTransViable; }
    public final boolean[] isTransCorrectLocation() { return mTransCorrectLocation; }
    public final int[] getExonsSkipped() { return mExonsSkipped; }

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
        if(mBreakends[FS_UP] == null && mBreakends[FS_DOWN] == null)
            return;

        if(mBreakends[FS_UP] != null && mBreakends[FS_DOWN] != null)
        {
            SvVarData varUp = mBreakends[FS_UP].getSV();
            SvCluster clusterUp = varUp.getCluster();
            SvVarData varDown = mBreakends[FS_DOWN].getSV();
            SvCluster clusterDown = varDown.getCluster();

            SvChain matchingChain = null;

            if(varUp != varDown && clusterUp == clusterDown)
            {
                // check for a matching chain if the clusters are the same
                matchingChain = clusterUp.findSameChainForSVs(varUp, varDown);

                if (matchingChain != null)
                {
                    final int chainData[] = breakendsAreChained(
                            matchingChain, varUp, !mBreakends[FS_UP].usesStart(), varDown, !mBreakends[FS_DOWN].usesStart());

                    mChainInfo = String.format("%d;%d", chainData[CHAIN_LINK_COUNT], chainData[CHAIN_LENGTH]);
                }
            }

            setFusionClusterInfo(mBreakends[FS_UP], true, matchingChain);
            setFusionClusterInfo(mBreakends[FS_DOWN], false, matchingChain);
        }
        else
        {
            SvBreakend breakend = mBreakends[FS_UP] != null ? mBreakends[FS_UP] : mBreakends[FS_DOWN];
            setFusionClusterInfo(breakend, mBreakends[FS_UP] != null, null);
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
