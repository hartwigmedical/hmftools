package com.hartwig.hmftools.common.purple.purity;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.stat.KernelDensity;
import org.jetbrains.annotations.NotNull;

class SomaticPeakFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticPeakFactory.class);
    private static final double KERNEL_BANDWIDTH = 0.02;

    public static List<SomaticPeak> findSomaticPeaks(@NotNull final List<? extends PurityAdjustedSomaticVariant> variants) {
        return findPeaks(variants.stream().map(PurityAdjustedSomaticVariant::adjustedVAF).collect(Collectors.toList()));
    }

    static List<SomaticPeak> findPeaks(@NotNull final List<Double> sample) {
        final SparkConf conf = new SparkConf().setAppName("com.hartwig.SomaticPeaks").setMaster("local");
        final JavaSparkContext sc = new JavaSparkContext(conf);
        final JavaRDD<Double> data = sc.parallelize(sample);
        final KernelDensity kd = new KernelDensity().setSample(data).setBandwidth(KERNEL_BANDWIDTH);

        final double[] vafs = IntStream.rangeClosed(0, 51).mapToDouble(x -> x / 100d).toArray();
        final double[] densities = kd.estimate(vafs);

        final List<SomaticPeak> results = Lists.newArrayList();
        for (int i = 1; i < densities.length - 1; i++) {
            double density = densities[i];
            if (Doubles.greaterThan(density, densities[i - 1]) && Doubles.greaterThan(density, densities[i + 1])) {
                final double adjustedVaf = vafs[i];
                final int peakCount = count(adjustedVaf, sample);
                final SomaticPeak peak = ImmutableSomaticPeak.builder().adjustedVAF(adjustedVaf).count(peakCount).build();
                LOGGER.info("Discovered peak {}", peak);
                results.add(peak);
            }
        }

        sc.stop();

        return results;
    }

    private static int count(double peak, @NotNull final List<Double> sample) {
        return (int) sample.stream().filter(vaf -> between(peak, vaf - KERNEL_BANDWIDTH / 2, vaf + KERNEL_BANDWIDTH / 2)).count();
    }

    private static boolean between(double victim, double min, double max) {
        return Doubles.greaterOrEqual(victim, min) && Doubles.lessOrEqual(victim, max);
    }
}
