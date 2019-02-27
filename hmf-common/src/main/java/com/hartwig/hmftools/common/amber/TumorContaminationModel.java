package com.hartwig.hmftools.common.amber;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class TumorContaminationModel {

    private static final Logger LOGGER = LogManager.getLogger(TumorContaminationModel.class);

    private static final double INCREMENT = 0.001;
    private static final long MIN_THREE_PLUS_READS = 2000;

    public double contamination(@NotNull final List<TumorContamination> samples) {
        long medianTumorReadDepth = medianDepth(samples);
        LOGGER.info("Median tumor depth at potential contamination sites is {} reads", medianTumorReadDepth);
        final Map<Integer, Long> map = samples.stream().collect(Collectors.groupingBy(x -> x.tumor().altSupport(), Collectors.counting()));
        return contamination(medianTumorReadDepth, map);
    }

    double contamination(final long medianTumorReadDepth, @NotNull final Map<Integer, Long> altSupportMap) {
        long three_plus_reads = reads(3, altSupportMap);
        if (three_plus_reads >= MIN_THREE_PLUS_READS) {

            double contamination = 0;
            double lowestScore = Double.MAX_VALUE;
            long two_plus_reads = reads(2, altSupportMap);

            for (double i = INCREMENT; Doubles.lessOrEqual(i, 1); i = i + INCREMENT) {
                double score = contaminationScore(i, two_plus_reads, medianTumorReadDepth, altSupportMap);
                if (Doubles.lessThan(score, lowestScore)) {
                    lowestScore = score;
                    contamination = i;
                }
            }

            LOGGER.warn("Found evidence of {}% contamination ", Math.floor(contamination * 1000) / 10);
            return contamination;
        }

        LOGGER.info("No evidence of contamination.");
        return 0;
    }

    private double contaminationScore(final double contamination, final long two_plus_reads, final long medianTumorReadDepth,
            @NotNull final Map<Integer, Long> altSupportMap) {

        final PoissonDistribution hetDistribution = new PoissonDistribution(0.5 * contamination * medianTumorReadDepth);
        final PoissonDistribution homAltDistribution = new PoissonDistribution(contamination * medianTumorReadDepth);

        final Function<Integer, Double> unadjustedModelLikelihood = altSupport -> {
            double hetLikelihood = hetDistribution.probability(altSupport);
            double homAltLikelihood = homAltDistribution.probability(altSupport);
            return 0.5 * hetLikelihood + 0.25 * homAltLikelihood;
        };

        final double altSupport_0 = unadjustedModelLikelihood.apply(0) + 0.25 * 1;
        final double altSupport_1 = unadjustedModelLikelihood.apply(1);
        final double altSupport_2_plus_adjustment = 1 - altSupport_0 - altSupport_1;

        final Function<Integer, Double> modelLikelihood =
                integer -> unadjustedModelLikelihood.apply(integer) / altSupport_2_plus_adjustment;

        double totalModelPercentage = 0;
        double totalDifference = 0;

        for (Map.Entry<Integer, Long> entry : altSupportMap.entrySet()) {
            final long altSupport = entry.getKey();
            final long altSupportCount = entry.getValue();

            if (altSupport > 1) {
                double modelPercentage = modelLikelihood.apply((int) altSupport);
                double actualPercentage = altSupportCount * 1d / two_plus_reads;

                totalDifference += Math.abs(actualPercentage - modelPercentage);
                totalModelPercentage += modelPercentage;
            }
        }

        return totalDifference + (1 - totalModelPercentage);
    }

    static long reads(int minAltSupport, @NotNull final Map<Integer, Long> altSupportMap) {
        return altSupportMap.entrySet().stream().filter(x -> x.getKey() >= minAltSupport).mapToLong(Map.Entry::getValue).sum();
    }

    private static int medianDepth(@NotNull final List<TumorContamination> baf) {
        final List<Integer> reads = baf.stream().map(x -> x.tumor().readDepth()).filter(x -> x > 0).sorted().collect(Collectors.toList());
        int count = reads.size();
        return count % 2 == 0 ? (reads.get(count / 2) + reads.get(count / 2 - 1)) / 2 : reads.get(count / 2);
    }

}
