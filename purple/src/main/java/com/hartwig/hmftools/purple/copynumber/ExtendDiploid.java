package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.function.IntUnaryOperator;

import com.google.common.base.MoreObjects;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public class ExtendDiploid
{
    private enum Direction
    {
        LEFT(index -> index - 1),
        RIGHT(index -> index + 1);

        private final IntUnaryOperator indexOperator;

        Direction(final IntUnaryOperator indexOperator)
        {
            this.indexOperator = indexOperator;
        }

        int moveIndex(int index)
        {
            return indexOperator.applyAsInt(index);
        }
    }

    private static final int MIN_BAF_COUNT_TO_WEIGH_WITH_BAF = 50;
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    private final int mMinTumorCount;
    private final int mCentromereMinTumorCount;
    private final CopyNumberTolerance mTolerance;

    public ExtendDiploid(final CopyNumberTolerance tolerance, final int minTumorCount, final int minTumorCountAtCentromere)
    {
        mMinTumorCount = minTumorCount;
        mCentromereMinTumorCount = minTumorCountAtCentromere;
        mTolerance = tolerance;
    }

    @NotNull
    public List<CombinedRegion> extendDiploid(final Collection<ObservedRegion> fittedRegions)
    {
        final boolean bafWeighted = fittedRegions.stream().anyMatch(x -> x.bafCount() >= MIN_BAF_COUNT_TO_WEIGH_WITH_BAF);

        final List<CombinedRegion> regions = Lists.newLinkedList();

        for(ObservedRegion fittedRegion : fittedRegions)
        {
            regions.add(new CombinedRegion(bafWeighted, fittedRegion));
        }

        int highestConfidenceIndex = nextIndex(regions);
        while(highestConfidenceIndex > -1)
        {
            final CombinedRegion highestConfidence = regions.get(highestConfidenceIndex);
            highestConfidence.setCopyNumberMethod(CopyNumberMethod.BAF_WEIGHTED);

            PPL_LOGGER.trace("Selected region {}", toString(highestConfidence.region()));
            extendRight(regions, highestConfidenceIndex);
            extendLeft(regions, highestConfidenceIndex);

            PPL_LOGGER.trace("Completed region {}", toString(highestConfidence.region()));
            highestConfidenceIndex = nextIndex(regions);
        }

        return regions;
    }

    private void extendRight(final List<CombinedRegion> regions, int targetIndex)
    {
        assert (targetIndex < regions.size());
        int neighbourIndex = targetIndex + 1;

        while(neighbourIndex < regions.size())
        {
            if(!merge(regions, Direction.RIGHT, targetIndex))
            {
                return;
            }
            regions.remove(neighbourIndex);
        }
    }

    private void extendLeft(final List<CombinedRegion> regions, final int targetIndex)
    {
        assert (targetIndex < regions.size());
        int neighbourIndex = targetIndex - 1;

        while(neighbourIndex >= 0)
        {
            if(!merge(regions, Direction.LEFT, neighbourIndex + 1))
            {
                return;
            }
            regions.remove(neighbourIndex);
            neighbourIndex--;
        }
    }

    private boolean merge(final List<CombinedRegion> regions, final Direction direction, int targetIndex)
    {
        final CombinedRegion target = regions.get(targetIndex);
        final ObservedRegion neighbour = regions.get(direction.moveIndex(targetIndex)).region();

        if(Extend.doNotExtend(target, neighbour))
            return false;

        int minTumorCount = nextBigBreakIsCentromereOrTelomere(regions, direction, targetIndex) ? mCentromereMinTumorCount : mMinTumorCount;

        final boolean isNeighbourDubious = isDubious(minTumorCount, neighbour);
        if(isNeighbourDubious)
        {
            if(inTolerance(target.region(), neighbour))
            {
                target.extendWithWeightedAverage(neighbour);
                return true;
            }
            else if(pushThroughDubiousRegion(minTumorCount, regions, direction, targetIndex))
            {
                target.extend(neighbour);
                return true;
            }
            else
            {
                return false;
            }
        }

        final boolean isNeighbourValid = isValid(minTumorCount, neighbour);
        if(!isNeighbourValid)
        {
            target.extend(neighbour);
            return true;
        }
        else if(inTolerance(target.region(), neighbour))
        {
            target.extendWithWeightedAverage(neighbour);
            return true;
        }

        return false;
    }

    private boolean isValid(int minTumorCount, final ObservedRegion region)
    {
        return region.germlineStatus() == GermlineStatus.DIPLOID && (region.support().isSV() || region.depthWindowCount() >= minTumorCount);
    }

    private boolean isDubious(int minTumorCount, final ObservedRegion region)
    {
        return region.germlineStatus() == GermlineStatus.DIPLOID && !region.support().isSV() && region.depthWindowCount() < minTumorCount;
    }

    private boolean nextBigBreakIsCentromereOrTelomere(
            final List<CombinedRegion> regions, final Direction direction, int targetIndex)
    {
        for(int i = direction.moveIndex(targetIndex); i >= 0 && i < regions.size(); i = direction.moveIndex(i))
        {
            final ObservedRegion neighbour = regions.get(i).region();

            if(neighbour.support() == SegmentSupport.CENTROMERE)
                return true;

            if(neighbour.support().isSV())
                return false;
        }

        return true;
    }

    private boolean pushThroughDubiousRegion(
            int minTumorCount, final List<CombinedRegion> regions, final Direction direction, int targetIndex)
    {
        int dubiousCount = 0;
        final CombinedRegion target = regions.get(targetIndex);
        for(int i = direction.moveIndex(targetIndex); i >= 0 && i < regions.size(); i = direction.moveIndex(i))
        {
            final ObservedRegion neighbour = regions.get(i).region();

            // Coming from left to right, EXCLUDE neighbour from decision on break.
            if(neighbour.start() > target.start())
            {
                if(neighbour.support() == SegmentSupport.CENTROMERE)
                {
                    return dubiousCount < minTumorCount;
                }
                if(neighbour.support().isSV())
                {
                    return false;
                }
            }

            boolean isDubious = isDubious(minTumorCount, neighbour);
            boolean inTolerance = inTolerance(target.region(), neighbour);

            if(isDubious && !inTolerance)
            {
                dubiousCount += neighbour.depthWindowCount();
                if(dubiousCount >= minTumorCount)
                {
                    return false;
                }
            }

            if(isValid(minTumorCount, neighbour))
            {
                return inTolerance;
            }

            // Coming from right to left, INCLUDE neighbour from decision on break.
            if(neighbour.start() < target.start())
            {
                if(neighbour.support() == SegmentSupport.CENTROMERE)
                {
                    return dubiousCount < minTumorCount;
                }
                if(neighbour.support().isSV())
                {
                    return false;
                }
            }
        }

        return dubiousCount < minTumorCount;
    }

    private boolean inTolerance(final ObservedRegion left, final ObservedRegion right)
    {
        return mTolerance.inTolerance(left, right);
    }

    private static int nextIndex(final List<CombinedRegion> regions)
    {
        int indexOfLargestBaf = -1;
        int indexOfLargestLength = -1;
        int indexOfTumorRatioCount = -1;

        int largestBAFCount = 0;
        long largestDepthWindowCount = 0;

        for(int i = 0; i < regions.size(); i++)
        {
            final CombinedRegion combined = regions.get(i);
            final ObservedRegion region = combined.region();
            if(!combined.isProcessed() && region.germlineStatus().equals(GermlineStatus.DIPLOID))
            {

                if(region.bafCount() > largestBAFCount)
                {
                    largestBAFCount = region.bafCount();
                    indexOfLargestBaf = i;
                }

                if(region.depthWindowCount() > largestDepthWindowCount)
                {
                    largestDepthWindowCount = region.depthWindowCount();
                    indexOfTumorRatioCount = i;
                }
            }
        }

        return indexOfLargestBaf > -1 ? indexOfLargestBaf : (indexOfTumorRatioCount > -1 ? indexOfTumorRatioCount : indexOfLargestLength);
    }

    private static String toString(final ObservedRegion region)
    {
        return MoreObjects.toStringHelper("FittedRegion")
                .omitNullValues()
                .add("chromosome", region.chromosome())
                .add("start", region.start())
                .add("end", region.end())
                .add("status", region.germlineStatus())
                .add("support", region.support())
                .add("copyNumber", FORMAT.format(region.tumorCopyNumber()))
                .toString();
    }
}
