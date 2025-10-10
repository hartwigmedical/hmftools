package com.hartwig.hmftools.purple.targeted;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionDataBuilder
{
    private final Map<String,List<TaggedRegion>> mTargetRegions;
    private final RefGenomeVersion mGenomeVersion;
    private final Map<String,List<PurpleSegment>> mPurpleSegments;
    private final Map<Chromosome,List<CobaltRatio>> mCobaltRatios;
    private final Map<String,List<PurpleCopyNumber>> mPurpleCopyNumbers;

    private final List<TargetRegionsCopyNumber> mTargetRegionData;

    public TargetRegionDataBuilder(
            final RefGenomeVersion genomeVersion, final Map<String,List<TaggedRegion>> targetRegions, final List<PurpleSegment> segments,
            final Map<Chromosome,List<CobaltRatio>> cobaltData, final List<PurpleCopyNumber> purpleCopyNumbers)
    {
        mGenomeVersion = genomeVersion;
        mCobaltRatios = cobaltData;

        mPurpleCopyNumbers = PurpleCopyNumber.buildChromosomeMap(purpleCopyNumbers);

        mTargetRegions = targetRegions;
        mPurpleSegments = PurpleSegment.buildChromosomeMap(segments);

        mTargetRegionData = Lists.newArrayList();
    }

    public List<TargetRegionsCopyNumber> targetRegionData() { return mTargetRegionData; }

    public void buildTargetRegionData()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            buildChromosomeData(chromosome);
        }
    }

    private void buildChromosomeData(final HumanChromosome chromosome)
    {
        String chrStr = mGenomeVersion.versionedChromosome(chromosome);

        List<TaggedRegion> targetRegions = mTargetRegions.get(chrStr);
        List<CobaltRatio> cobaltRatios = mCobaltRatios.get(chromosome);

        if(targetRegions == null || targetRegions.isEmpty() || cobaltRatios == null || cobaltRatios.isEmpty())
            return;

        List<PurpleCopyNumber> purpleCopyNumbers = mPurpleCopyNumbers.get(chrStr);
        List<PurpleSegment> purpleSegments = mPurpleSegments.get(chrStr);

        if(purpleCopyNumbers == null || purpleCopyNumbers.isEmpty() || purpleSegments == null || purpleSegments.isEmpty())
            return;

        PPL_LOGGER.trace("building chr({}) target region CN data from cobaltRatio({}) targetRegions({}) segments({}) copyNumbers({})",
                chrStr, cobaltRatios.size(), targetRegions.size(), purpleSegments.size(), purpleCopyNumbers.size());

        TaggedRegion currentTargetRegion = targetRegions.get(0);
        int targetRegionIndex = 0;

        PurpleSegment currentSegment = purpleSegments.get(0);
        int segmentIndex = 0;

        PurpleCopyNumber currentCopyNumber = purpleCopyNumbers.get(0);
        int copyNumberIndex = 0;

        for(CobaltRatio cobaltRatio : cobaltRatios)
        {
            ChrBaseRegion cobaltRegion = cobaltRatio.window();

            if(currentTargetRegion == null || !overlapsTargetRegion(cobaltRegion, currentTargetRegion))
            {
                int newIndex = findTargetRegionIndex(targetRegions, targetRegionIndex, cobaltRegion);

                if(newIndex == INVALID_DATA_INDEX)
                    continue;

                targetRegionIndex = newIndex;
                currentTargetRegion = targetRegions.get(targetRegionIndex);
            }

            if(currentCopyNumber == null || !overlapsPurpleCopyNumber(cobaltRegion, currentCopyNumber))
            {
                int newIndex = findPurpleCopyNumberIndex(purpleCopyNumbers, copyNumberIndex, cobaltRegion);

                if(newIndex == INVALID_DATA_INDEX)
                    continue;

                copyNumberIndex = newIndex;
                currentCopyNumber = purpleCopyNumbers.get(copyNumberIndex);
            }

            if(currentSegment == null || !overlapsPurpleSegment(cobaltRegion, currentSegment))
            {
                int newIndex = findPurpleSegmentIndex(purpleSegments, segmentIndex, cobaltRegion);

                if(newIndex == INVALID_DATA_INDEX)
                {
                    currentSegment = null; // can this be null or not?
                    // continue;
                }
                else
                {
                    segmentIndex = newIndex;
                    currentSegment = purpleSegments.get(segmentIndex);
                }
            }

            // now have all matching data, compile it
            List<TaggedRegion> overlappingTargetRegions = findOverlappingTargetRegions(targetRegions, targetRegionIndex, cobaltRegion);

            // could also check for a window overlapping multiple copy numbers but given the typical size difference this isn't worth it

            TargetRegionsCopyNumber targetRegionData = new TargetRegionsCopyNumber(
                    cobaltRatio, overlappingTargetRegions, currentCopyNumber,
                    currentSegment != null ? currentSegment.GermlineState : GermlineStatus.UNKNOWN);

            mTargetRegionData.add(targetRegionData);
        }
    }

    private static final int INVALID_DATA_INDEX = -1;

    private static boolean overlapsPurpleCopyNumber(final ChrBaseRegion cobaltRegion, final PurpleCopyNumber copyNumber)
    {
        return cobaltRegion.overlaps(copyNumber.chromosome(), copyNumber.start(), copyNumber.end());
    }

    private static int findPurpleCopyNumberIndex(
            final List<PurpleCopyNumber> copyNumbers, int currentIndex, final ChrBaseRegion cobaltRegion)
    {
        if(copyNumbers.get(currentIndex).start() > cobaltRegion.end())
            return INVALID_DATA_INDEX;

        ++currentIndex;
        while(currentIndex < copyNumbers.size())
        {
            PurpleCopyNumber copyNumber = copyNumbers.get(currentIndex);

            if(overlapsPurpleCopyNumber(cobaltRegion, copyNumber))
                return currentIndex;

            if(copyNumber.start() > cobaltRegion.end())
                break;

            ++currentIndex;
        }

        return INVALID_DATA_INDEX;
    }

    private static boolean overlapsTargetRegion(final ChrBaseRegion cobaltRegion, final TaggedRegion targetRegion)
    {
        return cobaltRegion.overlaps(targetRegion);
    }

    private static List<TaggedRegion> findOverlappingTargetRegions(
            final List<TaggedRegion> targetRegions, final int currentIndex, final ChrBaseRegion cobaltRegion)
    {
        List<TaggedRegion> regions = Lists.newArrayList(targetRegions.get(currentIndex));

        for(int i = 0; i <= 1; ++i)
        {
            boolean searchUp = (i == 0);
            int index = currentIndex + (searchUp ? 1 : -1);

            while(index >= 0 && index < targetRegions.size())
            {
                TaggedRegion region = targetRegions.get(index);

                if(overlapsTargetRegion(cobaltRegion, region))
                    regions.add(region);
                else
                    break;

                index += searchUp ? 1 : -1;
            }
        }

        return regions;
    }

    private static int findTargetRegionIndex(
            final List<TaggedRegion> targetRegions, int currentIndex, final ChrBaseRegion cobaltRegion)
    {
        if(targetRegions.get(currentIndex).start() > cobaltRegion.end())
            return INVALID_DATA_INDEX;

        ++currentIndex;
        while(currentIndex < targetRegions.size())
        {
            TaggedRegion region = targetRegions.get(currentIndex);

            if(overlapsTargetRegion(cobaltRegion, region))
                return currentIndex;

            if(region.start() > cobaltRegion.end())
                break;

            ++currentIndex;
        }

        return INVALID_DATA_INDEX;
    }

    private static boolean overlapsPurpleSegment(final ChrBaseRegion cobaltRegion, final PurpleSegment segment)
    {
        return cobaltRegion.overlaps(segment.chromosome(), segment.start(), segment.end());
    }

    private static int findPurpleSegmentIndex(
            final List<PurpleSegment> segments, int currentIndex, final ChrBaseRegion cobaltRegion)
    {
        if(segments.get(currentIndex).start() > cobaltRegion.end())
            return INVALID_DATA_INDEX;

        ++currentIndex;
        while(currentIndex < segments.size())
        {
            PurpleSegment segment = segments.get(currentIndex);

            if(overlapsPurpleSegment(cobaltRegion, segment))
                return currentIndex;

            if(segment.start() > cobaltRegion.end())
                break;

            ++currentIndex;
        }

        return INVALID_DATA_INDEX;
    }
}
