package com.hartwig.hmftools.svtools.fusion_likelihood;

import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.UNKNOWN;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.ChromosomeArm;

public class GeneRangeData
{
    public final com.hartwig.hmftools.common.gene.GeneData GeneData;
    public final ChromosomeArm Arm;
    public final String ChromosomeArm;

    // a set of merged and extended regions with matching phase and pre-gene combinations
    private List<GenePhaseRegion> mPhaseRegions;

    // intronic regions for each transcript
    private List<GenePhaseRegion> mTranscriptPhaseRegions;

    private Boolean mStreamUpOnly; // if a gene is forced to only be a 5' or 3' partner

    // maps from the DEL or DUP bucket length array index to overlap count
    private Map<Integer,Long> mDelFusionBaseCounts;
    private Map<Integer,Long> mDupFusionBaseCounts;

    private long[] mBaseOverlapCountDownstream;
    private long[] mBaseOverlapCountUpstream;

    public static final int NON_PROX_TYPE_SHORT_INV = 0;
    public static final int NON_PROX_TYPE_MEDIUM_INV = 1;
    public static final int NON_PROX_TYPE_LONG_SAME_ARM = 2;
    public static final int NON_PROX_TYPE_REMOTE = 3;
    public static final int NON_PROX_TYPE_MAX = NON_PROX_TYPE_REMOTE + 1;

    public GeneRangeData(final com.hartwig.hmftools.common.gene.GeneData geneData)
    {
        GeneData = geneData;
        mPhaseRegions = Lists.newArrayList();
        mTranscriptPhaseRegions = Lists.newArrayList();

        Arm = getChromosomalArm(geneData.Chromosome, geneData.GeneStart);
        ChromosomeArm = makeChrArmStr(geneData.Chromosome, Arm);

        mStreamUpOnly = null;

        mDelFusionBaseCounts = Maps.newHashMap();
        mDupFusionBaseCounts = Maps.newHashMap();

        mBaseOverlapCountUpstream = new long[NON_PROX_TYPE_MAX];
        mBaseOverlapCountDownstream = new long[NON_PROX_TYPE_MAX];
    }

    private static ChromosomeArm getChromosomalArm(final String chromosome, final int position)
    {
        final Integer centromerePos = RefGenomeCoordinates.COORDS_37.Centromeres.get(HumanChromosome.fromString(chromosome));

        if(centromerePos == null)
            return UNKNOWN;

        return position < centromerePos ? P_ARM : Q_ARM;
    }

    private static String makeChrArmStr(final String chr, final ChromosomeArm arm)
    {
        return chr + "_" + com.hartwig.hmftools.common.purple.ChromosomeArm.asStr(arm);
    }

    public final List<GenePhaseRegion> getPhaseRegions() { return mPhaseRegions; }
    public void setPhaseRegions(List<GenePhaseRegion> regions) { mPhaseRegions = regions; }

    public final List<GenePhaseRegion> getTranscriptPhaseRegions() { return mTranscriptPhaseRegions; }
    public void setTranscriptPhaseRegions(List<GenePhaseRegion> regions) { mTranscriptPhaseRegions = regions; }

    public boolean isStreamRestricted() { return mStreamUpOnly != null; }
    public boolean isOnlyUpstreamPartner() { return mStreamUpOnly != null ? mStreamUpOnly : false; }
    public boolean isOnlyDownstreamPartner() { return mStreamUpOnly != null ? !mStreamUpOnly : false; }
    public void setRestrictedStream(Boolean stream) { mStreamUpOnly = stream; }

    public Map<Integer,Long> getDelFusionBaseCounts() { return mDelFusionBaseCounts; }
    public Map<Integer,Long> getDupFusionBaseCounts() { return mDupFusionBaseCounts; }

    public boolean hasCodingTranscripts()
    {
        return mPhaseRegions.stream().anyMatch(x -> x.Phase != PHASE_NON_CODING);
    }

    public long getBaseOverlapCountUpstream(int type) { return mBaseOverlapCountUpstream[type]; }
    public void addBaseOverlapCountUpstream(int type, long count) { mBaseOverlapCountUpstream[type] += count; }
    public long getBaseOverlapCountDownstream(int type) { return mBaseOverlapCountDownstream[type]; }
    public void addBaseOverlapCountDownstream(int type, long count) { mBaseOverlapCountDownstream[type] += count; }

    public void clearOverlapCounts()
    {
        mDelFusionBaseCounts.clear();
        mDupFusionBaseCounts.clear();

        for(int i = 0; i < NON_PROX_TYPE_MAX; ++i)
        {
            mBaseOverlapCountUpstream[i] = 0;
            mBaseOverlapCountDownstream[i] = 0;
        }
    }

    public int phasedRegionTotal()
    {
        final boolean[] codingPhases = {false, false, true, true, true};

        int codingBases = (int)mPhaseRegions.stream()
                .filter(x -> x.hasAnyPhaseMatch(codingPhases))
                .mapToLong(x -> x.length())
                .sum();

        return codingBases;
    }

    public boolean hasProteinCoding()
    {
        return mPhaseRegions.stream().anyMatch(GenePhaseRegion::proteinCoding);
    }

    public String toString()
    {
        return String.format("%s:%s %s %d -> %d", GeneData.GeneId, GeneData.GeneName, ChromosomeArm, GeneData.GeneStart, GeneData.GeneEnd);
    }

}
