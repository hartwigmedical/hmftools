package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Deprecated
class PurpleCopyNumberSmoothing {

    @NotNull
    static List<PurpleCopyNumber> smooth(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        return smoothUnknownBAF(copyNumbers);
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
