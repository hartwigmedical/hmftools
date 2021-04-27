package com.hartwig.hmftools.purple.copynumber;

import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

class ExtendDiploidBAF
{
    // LOH Parameters
    private static final double LOH_COPY_NUMBER = 0.5;
    private static final long LOH_SMALL_REGION_SIZE = 1000;
    private static final long LOH_MAX_NEAREST_DISTANCE = 1_000_000;

    private static final long TINY_REGION_SIZE = 30;
    private static final double MIN_COPY_NUMBER_CHANGE = 0.5;
    private static final InferRegion INVALID_PAIR = new InferRegion(-1, -1, -1, -1);
    private static final Set<SegmentSupport> IGNORE_SUPPORT =
            EnumSet.of(SegmentSupport.CENTROMERE, SegmentSupport.TELOMERE, SegmentSupport.UNKNOWN);

    private final Map<Long, Long> mSimpleDupMap = Maps.newHashMap();

    ExtendDiploidBAF(@NotNull final List<StructuralVariant> simpleVariants)
    {
        for(StructuralVariant simpleVariant : simpleVariants)
        {
            StructuralVariantLeg end = simpleVariant.end();
            if(simpleVariant.type() == StructuralVariantType.DUP && end != null)
            {
                mSimpleDupMap.put(simpleVariant.start().cnaPosition(), end.cnaPosition());
            }
        }
    }

    @VisibleForTesting
    ExtendDiploidBAF(@NotNull final Map<Long, Long> simpleDupMap)
    {
        this.mSimpleDupMap.putAll(simpleDupMap);
    }

    @NotNull
    List<CombinedRegion> extendBAF(@NotNull final List<CombinedRegion> regions)
    {
        InferRegion inferRegion = nextRegion(false, regions);
        while(inferRegion.isValid())
        {
            inferBetween(inferRegion, regions);
            inferRegion = nextRegion(false, regions);
        }

        inferRegion = nextRegion(true, regions);
        while(inferRegion.isValid())
        {
            inferBetween(inferRegion, regions);
            inferRegion = nextRegion(true, regions);
        }

        return regions;
    }

    static double bafForTargetAllele(double targetAllele, double copyNumber)
    {
        if(Doubles.lessOrEqual(copyNumber, 1))
        {
            return 1;
        }

        double major = Math.min(targetAllele, copyNumber);
        double minor = copyNumber - major;

        return Math.max(major, minor) / (major + minor);
    }

    private static double singleSourceTargetPloidy(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        assert (inferRegion.isLeftValid() ^ inferRegion.isRightValid());

        final FittedRegion source;
        final Optional<FittedRegion> optionalDifferentSource;

        if(inferRegion.isLeftValid())
        {
            source = regions.get(inferRegion.leftSourceIndex).region();
            optionalDifferentSource = lookLeft(inferRegion, regions);
        }
        else
        {
            source = regions.get(inferRegion.rightSourceIndex).region();
            optionalDifferentSource = lookRight(inferRegion, regions);
        }

        if(isSingleTinyRegionWithGreaterCopyNumber(inferRegion, regions, source))
        {
            return source.minorAlleleCopyNumber();
        }

        if(optionalDifferentSource.isPresent())
        {
            final FittedRegion nearestDifferentSource = optionalDifferentSource.get();

            if(isMinorAlleleDifferent(source, nearestDifferentSource))
            {
                return source.majorAlleleCopyNumber();
            }
        }
        return source.minorAlleleCopyNumber();
    }

    private double multiSourceTargetPloidy(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        final FittedRegion primarySource;
        final FittedRegion secondarySource;
        if(regions.get(inferRegion.leftSourceIndex).bafCount() > regions.get(inferRegion.rightSourceIndex).bafCount())
        {
            primarySource = regions.get(inferRegion.leftSourceIndex).region();
            secondarySource = regions.get(inferRegion.rightSourceIndex).region();
        }
        else
        {
            primarySource = regions.get(inferRegion.rightSourceIndex).region();
            secondarySource = regions.get(inferRegion.leftSourceIndex).region();
        }

        if(minCopyNumberLessThanSourceMajorAllelePloidies(inferRegion, regions))
        {
            return primarySource.minorAlleleCopyNumber();
        }

        if(isMinorAlleleDifferent(primarySource, secondarySource) || isMajorAlleleDifferent(primarySource, secondarySource))
        {
            return minorOrMajorMovedTargetPloidy(primarySource, primarySource, secondarySource);
        }

        if(isSingleTinyRegionWithGreaterCopyNumber(inferRegion, regions, primarySource))
        {
            return primarySource.minorAlleleCopyNumber();
        }

        if(isSimpleDupSurroundedByLOH(inferRegion, regions))
        {
            return primarySource.minorAlleleCopyNumber();
        }

        final Optional<FittedRegion> optionalNearestDifferentSource = nearestDifferentRegion(inferRegion, regions);
        if(optionalNearestDifferentSource.isPresent())
        {
            final FittedRegion nearestSource = optionalNearestDifferentSource.get();

            if(isSingleSmallRegionFlankedByLargeLOH(nearestSource, inferRegion, regions))
            {
                return primarySource.minorAlleleCopyNumber();
            }

            final FittedRegion leftSource = regions.get(inferRegion.leftSourceIndex).region();
            final FittedRegion rightSource = regions.get(inferRegion.rightSourceIndex).region();

            // Need other to be the one we compared with when finding nearest different region
            final FittedRegion other = nearestSource.start() > rightSource.start() ? rightSource : leftSource;

            return minorOrMajorMovedTargetPloidy(primarySource, nearestSource, other);
        }

        return primarySource.minorAlleleCopyNumber();
    }

    @VisibleForTesting
    static double minorOrMajorMovedTargetPloidy(@NotNull final FittedRegion source, @NotNull final FittedRegion primary,
            @NotNull final FittedRegion secondary)
    {
        boolean isMinorDifferent = isMinorAlleleDifferent(primary, secondary);
        boolean isMajorDifferent = isMajorAlleleDifferent(primary, secondary);
        assert (isMinorDifferent || isMajorDifferent);

        if(isMinorDifferent && isMajorDifferent)
        {

            if(Doubles.lessThan(Math.abs(primary.minorAlleleCopyNumber() - secondary.majorAlleleCopyNumber()), MIN_COPY_NUMBER_CHANGE))
            {
                double average = (primary.minorAlleleCopyNumber() + secondary.majorAlleleCopyNumber()) / 2;
                return closestAllele(average, source);
            }

            if(Doubles.lessThan(Math.abs(primary.majorAlleleCopyNumber() - secondary.minorAlleleCopyNumber()), MIN_COPY_NUMBER_CHANGE))
            {
                double average = (primary.majorAlleleCopyNumber() + secondary.minorAlleleCopyNumber()) / 2;
                return closestAllele(average, source);
            }
        }

        if(isMajorDifferent)
        {
            // Correct
            return source.minorAlleleCopyNumber();
        }

        // Correct
        return source.majorAlleleCopyNumber();
    }

    private static double closestAllele(double target, @NotNull final FittedRegion source)
    {
        double sourceMinorDistanceFromAverage = Math.abs(source.minorAlleleCopyNumber() - target);
        double sourceMajorDistanceFromAverage = Math.abs(source.majorAlleleCopyNumber() - target);
        return Doubles.lessThan(sourceMinorDistanceFromAverage, sourceMajorDistanceFromAverage)
                ? source.minorAlleleCopyNumber()
                : source.majorAlleleCopyNumber();
    }

    private void inferBetween(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        assert (inferRegion.isValid());

        // Exactly one source available (XOR)
        final double targetPloidy = inferRegion.isLeftValid() ^ inferRegion.isRightValid()
                ? singleSourceTargetPloidy(inferRegion, regions)
                : multiSourceTargetPloidy(inferRegion, regions);

        extend(targetPloidy, inferRegion, regions);
    }

    @NotNull
    private static Optional<FittedRegion> lookRight(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        final FittedRegion rightSource = regions.get(inferRegion.rightSourceIndex).region();

        for(int i = inferRegion.rightSourceIndex + 1; i < regions.size(); i++)
        {
            final FittedRegion right = regions.get(i).region();
            if(right.support().equals(SegmentSupport.CENTROMERE))
            {
                return Optional.empty();
            }

            if(right.bafCount() > 0 && !IGNORE_SUPPORT.contains(right.support()) && isEitherAlleleDifferent(right, rightSource))
            {
                return Optional.of(right);
            }
        }

        return Optional.empty();
    }

    @NotNull
    private static Optional<FittedRegion> lookLeft(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        final FittedRegion leftSource = regions.get(inferRegion.leftSourceIndex).region();
        for(int i = inferRegion.leftSourceIndex - 1; i >= 0; i--)
        {

            if(regions.get(i + 1).region().support().equals(SegmentSupport.CENTROMERE))
            {
                return Optional.empty();
            }

            final FittedRegion left = regions.get(i).region();
            if(left.bafCount() > 0 && !IGNORE_SUPPORT.contains(left.support()) && isEitherAlleleDifferent(left, leftSource))
            {
                return Optional.of(left);
            }
        }

        return Optional.empty();
    }

    @NotNull
    private static Optional<FittedRegion> nearestDifferentRegion(@NotNull final InferRegion inferRegion,
            @NotNull final List<CombinedRegion> regions)
    {

        final FittedRegion leftSource = regions.get(inferRegion.leftSourceIndex).region();
        final FittedRegion rightSource = regions.get(inferRegion.rightSourceIndex).region();

        Optional<FittedRegion> optionalRight = lookRight(inferRegion, regions);
        Optional<FittedRegion> optionalLeft = lookLeft(inferRegion, regions);
        if(!optionalRight.isPresent())
        {
            return optionalLeft;
        }

        if(!optionalLeft.isPresent())
        {
            return optionalRight;
        }

        long leftDistance = leftSource.end() - optionalLeft.get().start();
        long rightDistance = optionalRight.get().start() - rightSource.start();

        return leftDistance < rightDistance ? optionalLeft : optionalRight;

    }

    private static boolean isEitherAlleleDifferent(@NotNull final FittedRegion left, @NotNull final FittedRegion right)
    {
        return isMajorAlleleDifferent(left, right) || isMinorAlleleDifferent(left, right);
    }

    private static boolean isMajorAlleleDifferent(@NotNull final FittedRegion left, @NotNull final FittedRegion right)
    {
        return Doubles.greaterThan(Math.abs(left.majorAlleleCopyNumber() - right.majorAlleleCopyNumber()), MIN_COPY_NUMBER_CHANGE);
    }

    private static boolean isMinorAlleleDifferent(@NotNull final FittedRegion left, @NotNull final FittedRegion right)
    {
        return Doubles.greaterThan(Math.abs(left.minorAlleleCopyNumber() - right.minorAlleleCopyNumber()), MIN_COPY_NUMBER_CHANGE);
    }

    private static void extend(double targetAllele, @NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        assert (inferRegion.isLeftValid());

        for(int i = inferRegion.leftTargetIndex; i <= inferRegion.rightTargetIndex; i++)
        {
            final CombinedRegion target = regions.get(i);

            target.setInferredTumorBAF(bafForTargetAllele(targetAllele, target.tumorCopyNumber()));
        }
    }

    @NotNull
    static InferRegion nextRegion(boolean crossCentromere, @NotNull final List<CombinedRegion> regions)
    {
        int leftSource = -1;
        boolean shouldInfer = false;
        int leftTarget = Integer.MAX_VALUE;

        for(int i = 0; i < regions.size(); i++)
        {
            CombinedRegion suspect = regions.get(i);
            if(suspect.support() == SegmentSupport.CENTROMERE && !crossCentromere)
            {
                if(shouldInfer && leftSource > -1)
                {
                    return new InferRegion(leftSource, leftTarget, i - 1, -1);
                }
                else
                {
                    leftSource = -1;
                    shouldInfer = false;
                    leftTarget = Integer.MAX_VALUE;
                }
            }

            if(suspect.bafCount() == 0 && !suspect.isInferredBAF())
            {
                shouldInfer = true;
                leftTarget = Math.min(leftTarget, i);
            }
            else if(shouldInfer)
            {
                return new InferRegion(leftSource, leftTarget, i - 1, i);
            }
            else
            {
                leftSource = i;
            }
        }

        return shouldInfer ? new InferRegion(leftSource, leftTarget, regions.size() - 1, -1) : INVALID_PAIR;
    }

    private static boolean isSingleTinyRegionWithGreaterCopyNumber(@NotNull final InferRegion inferRegion,
            @NotNull final List<CombinedRegion> regions, @NotNull final FittedRegion source)
    {
        return isSingleUnknownRegion(inferRegion) && isSmallRegion(TINY_REGION_SIZE, inferRegion, regions)
                && Doubles.greaterThan(regions.get(inferRegion.leftTargetIndex).tumorCopyNumber(), source.tumorCopyNumber());
    }

    private static boolean isSingleSmallRegionFlankedByLargeLOH(@NotNull final FittedRegion nearestSource,
            @NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        long distance = nearestSource.start() > regions.get(inferRegion.leftTargetIndex).end() ? nearestSource.start() - regions.get(
                inferRegion.leftTargetIndex).end() : regions.get(inferRegion.leftTargetIndex).start() - nearestSource.end();

        return isSmallRegionFlankedByLOH(inferRegion, regions) && distance > LOH_MAX_NEAREST_DISTANCE;
    }

    private static boolean isSmallRegionFlankedByLOH(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        return isSmallRegion(LOH_SMALL_REGION_SIZE, inferRegion, regions) && isFlankedByLOH(inferRegion, regions);
    }

    private static boolean isFlankedByLOH(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        return isLOH(regions.get(inferRegion.leftSourceIndex)) && isLOH(regions.get(inferRegion.rightSourceIndex));
    }

    private static boolean isLOH(@NotNull final CombinedRegion region)
    {
        return Doubles.lessOrEqual(region.region().minorAlleleCopyNumber(), LOH_COPY_NUMBER);
    }

    private static boolean isSmallRegion(long maxSize, @NotNull final InferRegion inferRegion,
            @NotNull final List<CombinedRegion> regions)
    {
        return regions.get(inferRegion.rightTargetIndex).end() - regions.get(inferRegion.leftTargetIndex).start() + 1 <= maxSize;
    }

    private static boolean isSingleUnknownRegion(@NotNull final InferRegion inferRegion)
    {
        return inferRegion.leftTargetIndex == inferRegion.rightTargetIndex;
    }

    private static boolean minCopyNumberLessThanSourceMajorAllelePloidies(@NotNull InferRegion inferRegion,
            @NotNull final List<CombinedRegion> regions)
    {
        assert (inferRegion.isRightValid());
        assert (inferRegion.isLeftValid());

        double minCopyNumber = minTargetCopyNumber(inferRegion, regions);
        return Doubles.greaterThan(regions.get(inferRegion.leftSourceIndex).region().majorAlleleCopyNumber() - minCopyNumber,
                MIN_COPY_NUMBER_CHANGE) && Doubles.greaterThan(
                regions.get(inferRegion.rightSourceIndex).region().majorAlleleCopyNumber() - minCopyNumber, MIN_COPY_NUMBER_CHANGE);

    }

    private static double minTargetCopyNumber(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        double result = Double.MAX_VALUE;

        for(int i = inferRegion.leftTargetIndex; i <= inferRegion.rightTargetIndex; i++)
        {
            result = Math.min(result, regions.get(i).tumorCopyNumber());
        }

        return result;
    }

    @VisibleForTesting
    boolean isSimpleDupSurroundedByLOH(@NotNull final InferRegion inferRegion, @NotNull final List<CombinedRegion> regions)
    {
        boolean isValidIndexes = inferRegion.leftTargetIndex == inferRegion.rightTargetIndex && inferRegion.leftSourceIndex != -1
                && inferRegion.rightSourceIndex != -1;
        if(!isValidIndexes)
        {
            return false;
        }

        FittedRegion target = regions.get(inferRegion.leftTargetIndex).region();
        FittedRegion right = regions.get(inferRegion.rightSourceIndex).region();
        boolean isStartAndEndDup = target.support() == SegmentSupport.DUP && right.support() == SegmentSupport.DUP;
        if(!isStartAndEndDup)
        {
            return false;
        }

        boolean isTheSameDup = mSimpleDupMap.containsKey(target.start()) && mSimpleDupMap.get(target.start()) == right.start();
        if(!isTheSameDup)
        {
            return false;
        }

        FittedRegion left = regions.get(inferRegion.leftSourceIndex).region();

        return Doubles.lessThan(left.minorAlleleCopyNumber(), 0.5) && Doubles.lessThan(right.minorAlleleCopyNumber(), 0.5);
    }

    static class InferRegion
    {
        final int leftSourceIndex;
        final int leftTargetIndex;
        final int rightTargetIndex;
        final int rightSourceIndex;

        InferRegion(final int leftSourceIndex, final int leftTargetIndex, final int rightTargetIndex, final int rightSourceIndex)
        {
            this.leftSourceIndex = leftSourceIndex;
            this.leftTargetIndex = leftTargetIndex;
            this.rightTargetIndex = rightTargetIndex;
            this.rightSourceIndex = rightSourceIndex;
        }

        boolean isValid()
        {
            return isLeftValid() || isRightValid();
        }

        boolean isLeftValid()
        {
            return leftSourceIndex != -1;
        }

        boolean isRightValid()
        {
            return rightSourceIndex != -1;
        }
    }
}
