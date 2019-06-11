package com.hartwig.hmftools.svanalysis.fusion;

import static com.hartwig.hmftools.svanalysis.fusion.GenePhaseRegion.PHASE_NON_CODING;
import static com.hartwig.hmftools.svanalysis.fusion.GenePhaseRegion.REGION_TYPE_NON_CODING;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;

public class GeneRangeData
{
    public final EnsemblGeneData GeneData;
    public final String Arm;
    private List<GenePhaseRegion> mPhaseRegions;
    private List<GenePhaseRegion> mCombinedPhaseRegions;

    // maps from the DEL or DUP bucket length array index to overlap count
    private Map<Integer,Long> mDelFusionBaseCounts;
    private Map<Integer,Long> mDupFusionBaseCounts;

    private long mTranslocationOverlapCount;
    private long mLocalArmOverlapCount;

    public GeneRangeData(final EnsemblGeneData geneData)
    {
        GeneData = geneData;
        mPhaseRegions = Lists.newArrayList();
        mCombinedPhaseRegions = Lists.newArrayList();

        Arm = SvUtilities.getChromosomalArm(geneData.Chromosome, geneData.GeneStart);

        mDelFusionBaseCounts = Maps.newHashMap();
        mDupFusionBaseCounts = Maps.newHashMap();

        mTranslocationOverlapCount = 0;
        mLocalArmOverlapCount = 0;
    }

    public final List<GenePhaseRegion> getPhaseRegions() { return mPhaseRegions; }
    public void addPhaseRegions(List<GenePhaseRegion> regions) { mPhaseRegions.addAll(regions); }

    public final List<GenePhaseRegion> getCombinedPhaseRegions() { return mCombinedPhaseRegions; }
    public void setCombinedPhaseRegions(List<GenePhaseRegion> regions) { mCombinedPhaseRegions = regions; }

    public Map<Integer,Long> getDelFusionBaseCounts() { return mDelFusionBaseCounts; }
    public Map<Integer,Long> getDupFusionBaseCounts() { return mDupFusionBaseCounts; }

    public boolean hasCodingTranscripts()
    {
        return mPhaseRegions.stream().anyMatch(x -> x.Phase != PHASE_NON_CODING);
    }

    public long getTranslocationOverlapCount() { return mTranslocationOverlapCount; }
    public void addTranslocationOverlapCount(long count) { mTranslocationOverlapCount += count; }

    public long getLocalArmOverlapCount() { return mLocalArmOverlapCount; }
    public void addLocalArmOverlapCount(long count) { mLocalArmOverlapCount += count; }


}
