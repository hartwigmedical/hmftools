package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.kde.KernelEstimator;
import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ClonalityCutoffKernel {

    private static final Logger LOGGER = LogManager.getLogger(ClonalityCutoffKernel.class);
    private static final double KERNEL_BANDWIDTH = 0.05;

    public static double clonalCutoff(@NotNull final List<PurityAdjustedSomaticVariant> variants) {
        return findClonalCutoff(variants.stream().map(PurityAdjustedSomaticVariant::ploidy).collect(Collectors.toList()));
    }

    private static double findClonalCutoff(@NotNull final List<Double> sample) {
        final KernelEstimator estimator = new KernelEstimator(0.001, KERNEL_BANDWIDTH);
        sample.forEach(x -> estimator.addValue(x, 1.0D));

        LOGGER.info("Examining ploidy troughs for clonality cutoff");
        final double[] ploidies = IntStream.rangeClosed(0, 100).mapToDouble(x -> x / 100d).toArray();
        final double[] densities = DoubleStream.of(ploidies).map(estimator::getProbability).toArray();

        double trough = 0;
        for (int i = 1; i < densities.length - 1; i++) {
            double density = densities[i];

            if (Doubles.lessThan(density, densities[i - 1]) && Doubles.lessThan(density, densities[i + 1])) {
                final double ploidy = ploidies[i];
                int troughCount = count(ploidy, sample);
                LOGGER.info("Discovered trough at ploidy {} with count {}", ploidy, troughCount);
                trough = ploidy;
            }
        }

        if (Doubles.isZero(trough)) {
            LOGGER.info("Unable to find clonality cutoff.");
        } else {
            LOGGER.info("Determined clonality cutoff at ploidy {}", trough);
        }

        return trough;
    }

    private static int count(double peak, @NotNull final List<Double> sample) {
        return (int) sample.stream().filter(vaf -> between(peak, vaf - 0.03, vaf + 0.03)).count();
    }

    private static boolean between(double victim, double min, double max) {
        return Doubles.greaterOrEqual(victim, min) && Doubles.lessOrEqual(victim, max);
    }
}
