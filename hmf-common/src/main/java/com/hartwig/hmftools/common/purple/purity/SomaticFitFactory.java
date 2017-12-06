package com.hartwig.hmftools.common.purple.purity;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.variant.PurpleSomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class SomaticFitFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticFitFactory.class);

    private final int minPeak;

    SomaticFitFactory(final int minPeak) {
        this.minPeak = minPeak;
    }

    @NotNull
    Optional<FittedPurity> fromSomatics(@NotNull final List<FittedPurity> candidates, @NotNull final List<PurpleSomaticVariant> variants) {
        double minPurity = candidates.stream().mapToDouble(FittedPurity::purity).min().orElse(0);
        double maxPurity = candidates.stream().mapToDouble(FittedPurity::purity).max().orElse(1);

        LOGGER.info("Looking for peak allelic frequency between purity {} and {}.", minPurity, maxPurity);
        final List<SomaticPeak> peaks = SomaticPeakFactory.findSomaticPeaks(variants);

        for (int i = peaks.size() - 1; i >= 0; i--) {
            SomaticPeak peak = peaks.get(i);
            double impliedPurity = peak.alleleFrequency() * 2;
            if (Doubles.greaterOrEqual(impliedPurity, minPurity) && Doubles.lessOrEqual(impliedPurity, maxPurity) && peak.count() > minPeak) {
                LOGGER.info("Somatic implied purity: {}", impliedPurity);
                return Optional.of(closest(impliedPurity, candidates));
            }
        }

        LOGGER.info("Unable to determine somatic implied purity.");
        return Optional.empty();
    }

    @NotNull
    private static FittedPurity closest(double impliedPurity, @NotNull final List<FittedPurity> candidates) {
        FittedPurity closest = candidates.get(0);
        double closestDistance = Math.abs(closest.purity() - impliedPurity);
        for (int i = 1; i < candidates.size(); i++) {
            FittedPurity candidate = candidates.get(i);
            double candidateDistance = Math.abs(candidate.purity() - impliedPurity);
            if (Doubles.lessThan(candidateDistance, closestDistance)) {
                closest = candidate;
                closestDistance = candidateDistance;
            }
        }

        return closest;
    }
}
