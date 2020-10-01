package com.hartwig.hmftools.common.purple.purity;

import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class SomaticFitFactory {

    private static final Logger LOGGER = LogManager.getLogger(SomaticFitFactory.class);

    private final int minPeak;
    private final int minSomatics;
    private final double minPurity;
    private final double maxPurity;

    SomaticFitFactory(final int minPeak, final int minSomatics, final double minPurity, final double maxPurity) {
        this.minPeak = minPeak;
        this.minSomatics = minSomatics;
        this.minPurity = minPurity;
        this.maxPurity = maxPurity;
    }

    @NotNull
    Optional<FittedPurity> fromSomatics(@NotNull final List<FittedPurity> allCandidates, @NotNull final List<SomaticVariant> variants) {

        if (variants.size() < minSomatics) {
            LOGGER.info("Insufficient somatic variants within average tumor depth to find somatic fit.");
            return Optional.empty();
        }

        LOGGER.info("Looking for peak somatic allelic frequencies");
        final List<SomaticPeak> peaks = SomaticPeakFactory.findSomaticPeaks(variants);

        // First try and get the largest implied purity where peak count > minPeak
        int maxPeak = 0;
        for (int i = peaks.size() - 1; i >= 0; i--) {
            SomaticPeak peak = peaks.get(i);
            double impliedPurity = peak.alleleFrequency() * 2;
            if (inPurityRange(impliedPurity)) {
                if (peak.count() >= minPeak) {
                    Optional<FittedPurity> diploid = diploid(impliedPurity, allCandidates);
                    if (diploid.isPresent()) {
                        LOGGER.info("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    } else {
                        LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
                maxPeak = Math.max(maxPeak, peak.count());
            }
        }

        // Failing that, get the implied purity with the largest peak
        if (maxPeak > 0) {
            for (int i = peaks.size() - 1; i >= 0; i--) {
                SomaticPeak peak = peaks.get(i);
                double impliedPurity = peak.alleleFrequency() * 2;
                if (peak.count() == maxPeak) {
                    Optional<FittedPurity> diploid = diploid(impliedPurity, allCandidates);
                    if (diploid.isPresent()) {
                        LOGGER.info("Somatic implied purity: {}", impliedPurity);
                        return diploid;
                    } else {
                        LOGGER.warn("Unable to find diploid solution for implied purity: {}", impliedPurity);
                    }
                }
            }
        }

        LOGGER.info("Unable to determine somatic implied purity.");
        return Optional.empty();
    }

    private boolean inPurityRange(double impliedPurity) {
        return Doubles.greaterOrEqual(impliedPurity, minPurity) && Doubles.lessOrEqual(impliedPurity, maxPurity);
    }

    private static Optional<FittedPurity> diploid(double purity, List<FittedPurity> diploidCandidates ) {
        return diploidCandidates.stream().filter(x -> Doubles.equal(x.purity(), purity)).findFirst();
    }
}
