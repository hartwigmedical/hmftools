package com.hartwig.hmftools.purple.fitting;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class SomaticKernelDensityPeaks
{

    private static final Logger LOGGER = LogManager.getLogger(SomaticKernelDensityPeaks.class);
    private static final double KERNEL_BANDWIDTH = 0.03;

    @NotNull
    static List<SomaticPeak> findSomaticPeaks(@NotNull final List<? extends AllelicDepth> variants) {
        return findPeaks(variants.stream().map(AllelicDepth::alleleFrequency).collect(Collectors.toList()));
    }

    @NotNull
    static List<SomaticPeak> findPeaks(@NotNull final List<Double> sample) {
        final KernelEstimator estimator = new KernelEstimator(0.001, KERNEL_BANDWIDTH);
        sample.forEach(x -> estimator.addValue(x, 1.0D));

        final double[] vafs = IntStream.rangeClosed(0, 51).mapToDouble(x -> x / 100d).toArray();
        final double[] densities = DoubleStream.of(vafs).map(estimator::getProbability).toArray();

        final List<SomaticPeak> results = Lists.newArrayList();
        for (int i = 1; i < densities.length - 1; i++) {
            double density = densities[i];
            if (Doubles.greaterThan(density, densities[i - 1]) && Doubles.greaterThan(density, densities[i + 1])) {
                final double alleleFrequency = vafs[i];
                final int peakCount = count(alleleFrequency, sample);
                final SomaticPeak peak = ImmutableSomaticPeak.builder().alleleFrequency(alleleFrequency).count(peakCount).build();
                LOGGER.debug("Discovered peak {}", peak);
                results.add(peak);
            }
        }

        return results;
    }

    private static int count(double peak, @NotNull final List<Double> sample) {
        return (int) sample.stream().filter(vaf -> between(peak, vaf - 0.015, vaf + 0.015)).count();
    }

    private static boolean between(double victim, double min, double max) {
        return Doubles.greaterOrEqual(victim, min) && Doubles.lessOrEqual(victim, max);
    }
}
