package com.hartwig.hmftools.common.purple.copynumber;

import java.util.EnumSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

class ExtendDiploidBAF {

    private static final double MIN_COPY_NUMBER_CHANGE = 0.5;
    private static final InferRegion INVALID_PAIR = new InferRegion(-1, -1, -1, -1);
    private static Set<SegmentSupport> IGNORE_SUPPORT =
            EnumSet.of(SegmentSupport.CENTROMERE, SegmentSupport.TELOMERE, SegmentSupport.UNKNOWN);

    @NotNull
    static List<CombinedRegion> extendBAF(@NotNull final List<CombinedRegion> regions) {

        InferRegion inferRegion = nextRegion(false, regions);
        while (inferRegion.isValid()) {
            inferBetween(inferRegion, regions);
            inferRegion = nextRegion(false, regions);
        }

        inferRegion = nextRegion(true, regions);
        while (inferRegion.isValid()) {
            inferBetween(inferRegion, regions);
            inferRegion = nextRegion(true, regions);
        }

        return regions;
    }

    static double bafForTargetAllele(double targetAllele, double copyNumber) {
        if (Doubles.isZero(copyNumber)) {
            return 0;
        }

        if (Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }

        double major = Math.min(targetAllele, copyNumber);
        double minor = copyNumber - major;

        return Math.max(major, minor) / (major + minor);
    }

    private static double singleSourceTargetPloidy(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions) {
        assert (inferRegion.isLeftValid() ^ inferRegion.isRightValid());

        final FittedRegion source;
        final Optional<FittedRegion> optionalDifferentSource;

        if (inferRegion.isLeftValid()) {
            source = regions.get(inferRegion.leftSource).region();
            optionalDifferentSource = lookLeft(inferRegion, regions);
        } else {
            source = regions.get(inferRegion.rightSource).region();
            optionalDifferentSource = lookRight(inferRegion, regions);
        }

        if (optionalDifferentSource.isPresent()) {
            final FittedRegion nearestDifferentSource = optionalDifferentSource.get();

            if (isMinorAlleleDifferent(source, nearestDifferentSource)) {
                return source.majorAllelePloidy();
            }
        }
        return source.minorAllelePloidy();
    }

    private static double multiSourceTargetPloidy(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions) {

        FittedRegion primarySource;
        FittedRegion secondarySource;
        if (regions.get(inferRegion.leftSource).bafCount() > regions.get(inferRegion.rightSource).bafCount()) {
            primarySource = regions.get(inferRegion.leftSource).region();
            secondarySource = regions.get(inferRegion.rightSource).region();
        } else {
            primarySource = regions.get(inferRegion.rightSource).region();
            secondarySource = regions.get(inferRegion.leftSource).region();
        }

        if (isMinorAlleleDifferent(primarySource, secondarySource) || isMajorAlleleDifferent(primarySource, secondarySource)) {
            return minorOrMajorMovedTargetPloidy(primarySource, secondarySource);
        }

        final Optional<FittedRegion> optionalNearestDifferentSource = nearestDifferentRegion(inferRegion, regions);
        if (optionalNearestDifferentSource.isPresent()) {

            final FittedRegion leftSource = regions.get(inferRegion.leftSource).region();
            final FittedRegion rightSource = regions.get(inferRegion.rightSource).region();

            secondarySource = optionalNearestDifferentSource.get();
            primarySource = secondarySource.start() < leftSource.start() ? rightSource : leftSource;

            return minorOrMajorMovedTargetPloidy(primarySource, secondarySource);
        }

        return primarySource.minorAllelePloidy();
    }

    private static double minorOrMajorMovedTargetPloidy(@NotNull final FittedRegion primarySource, @NotNull final FittedRegion secondarySource) {
        boolean isMinorDifferent = isMinorAlleleDifferent(primarySource, secondarySource);
        boolean isMajorDifferent = isMajorAlleleDifferent(primarySource, secondarySource);
        assert (isMinorDifferent || isMajorDifferent);

        if (isMinorDifferent && isMajorDifferent) {

            if (Doubles.lessThan(Math.abs(primarySource.minorAllelePloidy() - secondarySource.majorAllelePloidy()),
                    MIN_COPY_NUMBER_CHANGE)) {
                return primarySource.minorAllelePloidy();
            }

            if (Doubles.lessThan(Math.abs(primarySource.majorAllelePloidy() - secondarySource.minorAllelePloidy()),
                    MIN_COPY_NUMBER_CHANGE)) {
                return primarySource.majorAllelePloidy();
            }
        }

        if (isMajorDifferent) {
            return primarySource.minorAllelePloidy();
        }

        return primarySource.majorAllelePloidy();
    }


    private static void inferBetween(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions) {
        assert (inferRegion.isValid());

        // JOBA: Exactly one source available (XOR)
        final double targetPloidy = inferRegion.isLeftValid() ^ inferRegion.isRightValid()
                ? singleSourceTargetPloidy(inferRegion, regions)
                : multiSourceTargetPloidy(inferRegion, regions);

        extend(targetPloidy, inferRegion, regions);
    }

    @NotNull
    private static Optional<FittedRegion> lookRight(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions) {
        final FittedRegion rightSource = regions.get(inferRegion.rightSource).region();

        for (int i = inferRegion.rightSource + 1; i < regions.size(); i++) {
            FittedRegion right = regions.get(i).region();
            if (right.bafCount() > 0 && !IGNORE_SUPPORT.contains(right.support()) && isEitherAlleleDifferent(right, rightSource)) {
                return Optional.of(right);
            }
        }

        return Optional.empty();
    }

    @NotNull
    private static Optional<FittedRegion> lookLeft(@NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions) {
        final FittedRegion leftSource = regions.get(inferRegion.leftSource).region();

        for (int i = inferRegion.leftSource - 1; i >= 0; i--) {
            FittedRegion left = regions.get(i).region();
            if (left.bafCount() > 0 && !IGNORE_SUPPORT.contains(left.support()) && isEitherAlleleDifferent(left, leftSource)) {
                return Optional.of(left);
            }
        }

        return Optional.empty();
    }

    @NotNull
    private static Optional<FittedRegion> nearestDifferentRegion(@NotNull InferRegion inferRegion,
            @NotNull final List<CombinedRegion> regions) {

        final FittedRegion leftSource = regions.get(inferRegion.leftSource).region();
        final FittedRegion rightSource = regions.get(inferRegion.rightSource).region();

        Optional<FittedRegion> optionalRight = lookRight(inferRegion, regions);
        Optional<FittedRegion> optionalLeft = lookLeft(inferRegion, regions);

        if (!optionalRight.isPresent()) {
            return optionalLeft;
        }

        if (!optionalLeft.isPresent()) {
            return optionalRight;
        }

        long leftDistance = leftSource.end() - optionalLeft.get().start();
        long rightDistance = optionalRight.get().start() - rightSource.start();

        return leftDistance < rightDistance ? optionalLeft : optionalRight;

    }

    private static boolean isEitherAlleleDifferent(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        return isMajorAlleleDifferent(left, right) || isMinorAlleleDifferent(left, right);
    }

    private static boolean isMajorAlleleDifferent(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        return Doubles.greaterThan(Math.abs(left.majorAllelePloidy() - right.majorAllelePloidy()), MIN_COPY_NUMBER_CHANGE);
    }

    private static boolean isMinorAlleleDifferent(@NotNull final FittedRegion left, @NotNull final FittedRegion right) {
        return Doubles.greaterThan(Math.abs(left.minorAllelePloidy() - right.minorAllelePloidy()), MIN_COPY_NUMBER_CHANGE);
    }

    private static void extend(double targetAllele, @NotNull InferRegion inferRegion, @NotNull final List<CombinedRegion> regions) {
        assert (inferRegion.isLeftValid());

        for (int i = inferRegion.leftTarget; i <= inferRegion.rightTarget; i++) {
            final CombinedRegion target = regions.get(i);

            target.setInferredTumorBAF(bafForTargetAllele(targetAllele, target.tumorCopyNumber()));
        }
    }

    @NotNull
    static InferRegion nextRegion(boolean crossCentromere, @NotNull final List<CombinedRegion> regions) {

        int leftScope = -1;
        boolean shouldInfer = false;
        int leftTarget = Integer.MAX_VALUE;

        for (int i = 0; i < regions.size(); i++) {
            CombinedRegion suspect = regions.get(i);
            if (suspect.support() == SegmentSupport.CENTROMERE && !crossCentromere) {
                if (shouldInfer && leftScope > -1) {
                    return new InferRegion(leftScope, leftTarget, i - 1, -1);
                } else {
                    leftScope = -1;
                    shouldInfer = false;
                    leftTarget = Integer.MAX_VALUE;
                }
            }

            if (suspect.bafCount() == 0 && !suspect.isInferredBAF()) {
                shouldInfer = true;
                leftTarget = Math.min(leftTarget, i);
            } else if (shouldInfer) {
                return new InferRegion(leftScope, leftTarget, i - 1, i);
            } else {
                leftScope = i;
            }
        }

        return shouldInfer ? new InferRegion(leftScope, leftTarget, regions.size() - 1, -1) : INVALID_PAIR;
    }

    static class InferRegion {
        final int leftSource;
        final int leftTarget;
        final int rightTarget;
        final int rightSource;

        InferRegion(final int leftSource, final int leftTarget, final int rightTarget, final int rightSource) {
            this.leftSource = leftSource;
            this.leftTarget = leftTarget;
            this.rightTarget = rightTarget;
            this.rightSource = rightSource;
        }

        boolean isValid() {
            return isLeftValid() || isRightValid();
        }

        boolean isLeftValid() {
            return leftSource != -1;
        }

        boolean isRightValid() {
            return rightSource != -1;
        }
    }
}
