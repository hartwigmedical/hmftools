package com.hartwig.hmftools.purple.fitting;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.Doubles.lessOrEqual;

import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.purple.purity.FittedPurityScoreFactory;
import com.hartwig.hmftools.common.purple.purity.ImmutableBestFit;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.purple.config.ConfigSupplier;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class BestFitFactory {

    private static final DecimalFormat FORMAT = new DecimalFormat("0%");
    private static final Logger LOGGER = LogManager.getLogger(BestFitFactory.class);

    private static final double PERCENT_RANGE = 0.1;
    private static final double ABS_RANGE = 0.0005;

    private final ConfigSupplier config;

    private final int minReadCount;
    private final int maxReadCount;
    private final boolean somaticFitEnabled;

    private final double highlyDiploidPercentage;
    private final double minSomaticPurity;
    private final double minSomaticPuritySpread;

    private final int minTotalSvFragmentCount;
    private final int minTotalSomaticVariantAlleleReadCount;

    private final SomaticPurityFitter somaticPurityFitter;

    @NotNull
    private final BestFit bestFit;

    public BestFitFactory(final ConfigSupplier configSupplier,
            boolean somaticFitEnabled, int minReadCount, int maxReadCount, double minPurity, double maxPurity,
            int minVariants, int minPeak, double highlyDiploidPercentage, double minSomaticPurity, double minSomaticPuritySpread,
            int minTotalSvFragmentCount, int minTotalSomaticVariantAlleleReadCount, @NotNull final List<FittedPurity> allCandidates,
            @NotNull final List<SomaticVariant> somatics, @NotNull final List<StructuralVariant> structuralVariants) {
        assert (!allCandidates.isEmpty());
        this.config = configSupplier;
        this.somaticFitEnabled = somaticFitEnabled;
        this.minReadCount = minReadCount;
        this.maxReadCount = maxReadCount;
        this.minSomaticPuritySpread = minSomaticPuritySpread;
        this.minTotalSvFragmentCount = minTotalSvFragmentCount;
        this.minTotalSomaticVariantAlleleReadCount = minTotalSomaticVariantAlleleReadCount;

        this.minSomaticPurity = minSomaticPurity;
        this.highlyDiploidPercentage = highlyDiploidPercentage;

        this.somaticPurityFitter = new SomaticPurityFitter(configSupplier, minPeak, minVariants, minPurity, maxPurity);
        this.bestFit = bestFit(allCandidates, somatics, structuralVariants);
    }

    private BestFit bestFit(@NotNull final List<FittedPurity> allCandidates, @NotNull final List<SomaticVariant> somatics,
            @NotNull final List<StructuralVariant> structuralVariants) {
        Collections.sort(allCandidates);
        FittedPurity lowestScoreFit = allCandidates.get(0);

        final List<FittedPurity> bestFitPerPurityCandidates = inRangeOfLowest(lowestScoreFit.score(), allCandidates);
        final FittedPurityScore score = FittedPurityScoreFactory.score(bestFitPerPurityCandidates);

        final ImmutableBestFit.Builder builder = ImmutableBestFit.builder().score(score).allFits(allCandidates);

        boolean useSomatics;

        if(config.commonConfig().tumorOnly())
        {
            useSomatics = false;
        }
        else if(config.somaticConfig().forceSomaticFit())
        {
            useSomatics = true;

            LOGGER.info("forcing somatic fit");
        }
        else
        {
            useSomatics = Doubles.greaterOrEqual(score.puritySpread(), minSomaticPuritySpread) && isHighlyDiploid(score);

            if(useSomatics)
            {
                LOGGER.info("Sample is highly diploid [{}] with large purity range [{}:{}]",
                        FORMAT.format(score.maxDiploidProportion()), FORMAT.format(score.minPurity()), FORMAT.format(score.maxPurity()));
            }
        }

        if (!useSomatics)
            return builder.fit(lowestScoreFit).method(FittedPurityMethod.NORMAL).build();

        final List<FittedPurity> diploidCandidates = BestFit.mostDiploidPerPurity(allCandidates);

        if (diploidCandidates.isEmpty())
        {
            LOGGER.warn("Unable to use somatic fit as there are no diploid candidates");
            return builder.fit(lowestScoreFit).method(FittedPurityMethod.NORMAL).build();
        }

        final FittedPurity lowestPurityFit = diploidCandidates.stream().min(Comparator.comparingDouble(FittedPurity::purity)).get();

        final SvSummary svSummary = new SvSummary(structuralVariants);
        final SomaticSummary somaticSummary = new SomaticSummary(somatics);
        boolean hasTumor = somaticSummary.hotspots() > 0 || somaticSummary.totalAlleleReadCount() >= minTotalSomaticVariantAlleleReadCount
                || svSummary.hotspots() > 0 || svSummary.totalFragmentCount() >= minTotalSvFragmentCount;

        if (!hasTumor) {
            LOGGER.warn("No tumor found");
            return builder.fit(lowestPurityFit).method(FittedPurityMethod.NO_TUMOR).build();
        }

        final Optional<FittedPurity> somaticFit = somaticPurityFitter.fromSomatics(diploidCandidates, somaticSummary.filteredVariants());
        if (!somaticFit.isPresent()) {
            return builder.fit(lowestPurityFit).method(FittedPurityMethod.SOMATIC).build();
        } else if (somaticFitIsWorse(lowestScoreFit, somaticFit.get())) {
            return builder.fit(lowestScoreFit).method(FittedPurityMethod.NORMAL).build();
        } else {
            return builder.fit(somaticFit.get()).method(FittedPurityMethod.SOMATIC).build();
        }
    }

    private boolean somaticFitIsWorse(@NotNull final FittedPurity lowestScore, @NotNull final FittedPurity somaticFit) {
        double lowestPurity = lowestScore.purity();
        double somaticPurity = somaticFit.purity();

        return Doubles.lessThan(lowestPurity, minSomaticPurity) && Doubles.lessThan(somaticPurity, minSomaticPurity) && Doubles.greaterThan(
                somaticPurity,
                lowestPurity);
    }

    private boolean isHighlyDiploid(@NotNull final FittedPurityScore score) {
        return Doubles.greaterOrEqual(score.maxDiploidProportion(), highlyDiploidPercentage);
    }

    @NotNull
    public BestFit bestFit() {
        return bestFit;
    }

    @NotNull
    private static List<FittedPurity> inRangeOfLowest(double lowestScore, @NotNull final List<FittedPurity> purities) {
        return purities.stream().filter(inRangeOfLowest(lowestScore)).collect(toList());
    }

    @NotNull
    private static Predicate<FittedPurity> inRangeOfLowest(final double score) {
        return fittedPurity -> {
            double absDifference = Math.abs(fittedPurity.score() - score);
            double relDifference = Math.abs(absDifference / score);
            return lessOrEqual(absDifference, ABS_RANGE) || lessOrEqual(relDifference, PERCENT_RANGE);
        };
    }

    static class SvSummary {

        private int hotspotCount = 0;
        private int fragmentReadCount = 0;

        public SvSummary(@NotNull final List<StructuralVariant> variants) {
            for (StructuralVariant variant : variants) {
                if (variant.isFiltered()) {
                    continue;
                }

                if (variant.hotspot()) {
                    hotspotCount++;
                }

                Integer startTumorVariantFragmentCount = variant.start().tumorVariantFragmentCount();
                if (variant.end() != null && startTumorVariantFragmentCount != null) {
                    fragmentReadCount += startTumorVariantFragmentCount;
                }
            }
        }

        public int hotspots() {
            return hotspotCount;
        }

        public int totalFragmentCount() {
            return fragmentReadCount;
        }

    }

    class SomaticSummary {

        private int hotspotCount = 0;
        private int alleleReadCountTotal = 0;
        private final List<SomaticVariant> variantsInReadCountRange = Lists.newArrayList();

        public SomaticSummary(@NotNull final List<SomaticVariant> somatics) {
            for (SomaticVariant somatic : somatics) {
                if (somatic.isFiltered() || !somatic.isSnp()) {
                    continue;
                }

                if (somatic.isHotspot()) {
                    hotspotCount++;
                }

                alleleReadCountTotal += somatic.alleleReadCount();

                if (somatic.totalReadCount() >= minReadCount && somatic.totalReadCount() <= maxReadCount) {
                    variantsInReadCountRange.add(somatic);
                }
            }
        }

        public int hotspots() {
            return hotspotCount;
        }

        public int totalAlleleReadCount() {
            return alleleReadCountTotal;
        }

        @NotNull
        public List<SomaticVariant> filteredVariants() {
            return variantsInReadCountRange;
        }
    }
}
