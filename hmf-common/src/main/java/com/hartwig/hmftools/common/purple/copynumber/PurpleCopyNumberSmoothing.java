package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class PurpleCopyNumberSmoothing {

    @NotNull
    static List<FittedRegion> smoothRegions(@NotNull final List<FittedRegion> regions) {
        if (regions.size() > 1) {
            final FittedRegion telomere = regions.get(0);
            final FittedRegion neighbour = regions.get(1);
            if (Doubles.isZero(telomere.tumorCopyNumber())) {
                final List<FittedRegion> result = Lists.newArrayList();
                final FittedRegion fixedTelomere =
                        ImmutableFittedRegion.builder().from(telomere).tumorCopyNumber(neighbour.tumorCopyNumber()).build();
                result.add(fixedTelomere);
                result.addAll(regions.subList(1, regions.size()));
                return result;
            }
        }

        return regions;
    }

    @NotNull
    static List<PurpleCopyNumber> smooth(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        return smoothUnknownBAF(smoothTelomeres(copyNumbers));
    }

    @NotNull
    static List<PurpleCopyNumber> smoothTelomeres(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        if (copyNumbers.size() > 1) {
            final PurpleCopyNumber telomere = copyNumbers.get(0);
            final PurpleCopyNumber neighbour = copyNumbers.get(1);
            if (Doubles.isZero(telomere.averageTumorCopyNumber())) {
                final PurpleCopyNumber fixed = ImmutablePurpleCopyNumber.builder()
                        .from(telomere)
                        .averageTumorCopyNumber(neighbour.averageTumorCopyNumber())
                        .averageActualBAF(neighbour.averageActualBAF())
                        .build();
                copyNumbers.set(0, fixed);
            }
        }

        return copyNumbers;
    }

    @NotNull
    static List<PurpleCopyNumber> smoothUnknownBAF(@NotNull final List<PurpleCopyNumber> copyNumbers) {

        for (int i = 0; i < copyNumbers.size(); i++) {

            final PurpleCopyNumber region = copyNumbers.get(i);
            if (!Doubles.isZero(region.averageTumorCopyNumber()) && Doubles.isZero(region.averageActualBAF())) {
                final PurpleCopyNumber prev = i > 0 ? copyNumbers.get(i - 1) : null;
                final PurpleCopyNumber next = i < copyNumbers.size() - 1 ? copyNumbers.get(i + 1) : null;
                if (prev != null || next != null) {
                    final PurpleCopyNumber closest = closest(region.averageTumorCopyNumber(), prev, next);
                    final PurpleCopyNumber fixed = ImmutablePurpleCopyNumber.builder()
                            .from(region)
                            .averageActualBAF(SmoothBAF.estimateBAF(region, closest))
                            .build();

                    copyNumbers.set(i, fixed);
                }
            }
        }

        return copyNumbers;
    }

    @NotNull
    static PurpleCopyNumber closest(double copyNumber, @Nullable PurpleCopyNumber previous, @Nullable PurpleCopyNumber next) {
        assert (previous != null || next != null);

        if (previous == null) {
            return next;
        } else if (next == null) {
            return previous;
        }

        double previousDifference = Math.abs(copyNumber - previous.averageTumorCopyNumber());
        double nextDifference = Math.abs(copyNumber - next.averageTumorCopyNumber());

        return Doubles.lessThan(previousDifference, nextDifference) ? previous : next;
    }

}
