package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionsCopyNumbers
{
    private final List<TargetRegionsCopyNumber> mCopyNumbers = new ArrayList<>();

    public TargetRegionsCopyNumbers(final TargetRegionsData mTargetRegionsData,
            final Map<Chromosome, List<CobaltRatio>> mCobaltData,
            final List<PurpleCopyNumber> purpleCopyNumbers,
            final RefGenomeVersion mRefGenomeVersion)
    {
        SortedMap<CobaltRatio, List<TaggedRegion>> relevantCobaltRegions = new TreeMap<>();
        mCobaltData.forEach((chromosome, cobaltRatios) ->
        {
            List<TaggedRegion> targetRegions = mTargetRegionsData.targetRegions(mRefGenomeVersion.versionedChromosome(chromosome));
            cobaltRatios.forEach(cobaltRatio ->
            {
                List<TaggedRegion> overlappingRegions = cobaltRatio.findWindowOverlaps(targetRegions);
                if(!overlappingRegions.isEmpty())
                {
                    Collections.sort(overlappingRegions);
                    relevantCobaltRegions.put(cobaltRatio, overlappingRegions);
                }
            });
        });

        relevantCobaltRegions.forEach((cobaltRatio, regions) ->
        {
            List<PurpleCopyNumber> purpleRegionsContainingCobaltPoint = cobaltRatio.findContainingRegions(purpleCopyNumbers);
            if(purpleRegionsContainingCobaltPoint.size() == 1)
            {
                TargetRegionsCopyNumber regionsCopyNumber = new TargetRegionsCopyNumber(cobaltRatio, regions, purpleRegionsContainingCobaltPoint.get(0));
                mCopyNumbers.add(regionsCopyNumber);
            }
            else
            {
                PPL_LOGGER.info("Zero or multiple purple regions for position {} on chromosome {}", cobaltRatio.position(), cobaltRatio.chromosome());
            }
        });
    }

    public List<TargetRegionsCopyNumber> copyNumbersData()
    {
        return mCopyNumbers;
    }
}
