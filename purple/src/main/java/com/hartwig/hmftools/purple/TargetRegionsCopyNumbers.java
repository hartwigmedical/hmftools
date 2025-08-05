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
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

import org.apache.commons.lang3.tuple.Pair;

public class TargetRegionsCopyNumbers
{
    private final List<TargetRegionsCopyNumber> mCopyNumbers = new ArrayList<>();

    public TargetRegionsCopyNumbers(
            final TargetRegionsData mTargetRegionsData,
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
            List<PurpleCopyNumber> overlaps = cobaltRatio.window().findOverlaps(purpleCopyNumbers);
            if(overlaps.isEmpty())
            {
                PPL_LOGGER.info("No Purple regions for position {} on chromosome {}", cobaltRatio.position(), cobaltRatio.chromosome());
            }
            else
            {
                List<Pair<ChrBaseRegion, PurpleCopyNumber>> coverMap = cobaltRatio.window().splitByOverlappingRegions(overlaps);
                for(Pair<ChrBaseRegion, PurpleCopyNumber> pair : coverMap)
                {
                    final ChrBaseRegion region = pair.getLeft();
                    CobaltRatio restrictedRatio = cobaltRatio.realign(region.start());
                    final PurpleCopyNumber mPurpleCopyNumber = pair.getValue();
                    TargetRegionsCopyNumber regionsCopyNumber = new TargetRegionsCopyNumber(restrictedRatio, regions, mPurpleCopyNumber);
                    mCopyNumbers.add(regionsCopyNumber);
                }
            }
        });
    }

    public List<TargetRegionsCopyNumber> copyNumbersData()
    {
        return mCopyNumbers;
    }
}
