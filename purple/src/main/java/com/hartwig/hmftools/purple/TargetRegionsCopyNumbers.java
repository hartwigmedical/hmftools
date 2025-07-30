package com.hartwig.hmftools.purple;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionsCopyNumbers
{
    private final TargetRegionsData mTargetRegionsData;
    private final Map<Chromosome, List<CobaltRatio>> mCobaltData;
    private final Map<CobaltRatio, List<TaggedRegion>> mCobaltRegions = new HashMap<>();
    private final RefGenomeVersion mRefGenomeVersion;

    public TargetRegionsCopyNumbers(final TargetRegionsData mTargetRegionsData, final Map<Chromosome, List<CobaltRatio>> mCobaltData,
            final RefGenomeVersion mRefGenomeVersion)
    {
        this.mTargetRegionsData = mTargetRegionsData;
        this.mCobaltData = mCobaltData;
        this.mRefGenomeVersion = mRefGenomeVersion;
        findCobaltOverlaps();
    }

    private void findCobaltOverlaps()
    {
        mCobaltData.forEach((chromosome, cobaltRatios) ->
        {
            List<TaggedRegion> targetRegions = mTargetRegionsData.targetRegions(mRefGenomeVersion.versionedChromosome(chromosome));

            cobaltRatios.forEach(cobaltRatio ->
            {
                final List<TaggedRegion> overlappingRegions = new ArrayList<>();
                targetRegions.forEach(region ->
                {
                    if(cobaltRatio.baseRegion().overlaps(region))
                    {
                        overlappingRegions.add(region);
                    }
                });
                if(!overlappingRegions.isEmpty())
                {
                    mCobaltRegions.put(cobaltRatio, overlappingRegions);
                }
            });
        });
    }

    public List<TargetRegionsCopyNumber> copyNumbersData()
    {
        List<TargetRegionsCopyNumber> copyNumbers = new ArrayList<>();
        mCobaltRegions.forEach((cobaltRatio, regions) ->
                copyNumbers.add(new TargetRegionsCopyNumber(cobaltRatio, regions)));
        return copyNumbers;
    }
}
