package com.hartwig.hmftools.purple.region;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurityAdjusterAbnormalChromosome;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ImmutableFittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class FittedRegionFactoryV2 implements FittedRegionFactory {

    private final double ambiguousBaf;
    private final double ploidyPenaltyFactor;
    private final PloidyDeviation ploidyDeviation;
    private final CobaltChromosomes cobaltChromosomes;

    public FittedRegionFactoryV2(final CobaltChromosomes cobaltChromosomes, final int averageReadDepth, double ploidyPenaltyFactor,
            double ploidyPenaltyStandardDeviation, double ploidyPenaltyMinStandardDeviationPerPloidy,
            final double majorAlleleSubOnePenaltyMultiplier, final double majorAlleleSubOneAdditionalPenalty,
            final double baselineDeviation) {
        this.cobaltChromosomes = cobaltChromosomes;
        this.ploidyPenaltyFactor = ploidyPenaltyFactor;
        ploidyDeviation = new PloidyDeviation(ploidyPenaltyStandardDeviation,
                ploidyPenaltyMinStandardDeviationPerPloidy,
                majorAlleleSubOnePenaltyMultiplier,
                majorAlleleSubOneAdditionalPenalty,
                baselineDeviation);
        ambiguousBaf = ExpectedBAF.expectedBAF(averageReadDepth);
    }

    @Override
    @NotNull
    public List<FittedRegion> fitRegion(double purity, double normFactor, @NotNull final Collection<ObservedRegion> observedRegions) {
        final Predicate<ObservedRegion> valid = observedRegion -> isAllowedRegion(cobaltChromosomes, observedRegion);
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean isAllowedRegion(@NotNull final CobaltChromosomes cobaltChromosomes, @NotNull final GenomeRegion region) {
        return cobaltChromosomes.contains(region.chromosome());
    }

    @Override
    @NotNull
    public FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion) {
        final PurityAdjuster purityAdjuster = new PurityAdjusterAbnormalChromosome(purity, normFactor, cobaltChromosomes.chromosomes());

        double observedTumorRatio = observedRegion.observedTumorRatio();
        double impliedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double observedBAF = observedRegion.observedBAF();
        double impliedBAF = impliedBaf(purityAdjuster, observedRegion.chromosome(), impliedCopyNumber, observedBAF);

        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());

        double majorAllelePloidy = impliedBAF * impliedCopyNumber;
        double minorAllelePloidy = impliedCopyNumber - majorAllelePloidy;

        double majorAllelePloidyDeviation = ploidyDeviation.majorAlleleDeviation(purity, normFactor, majorAllelePloidy);
        double minorAllelePloidyDeviation = ploidyDeviation.minorAlleleDeviation(purity, normFactor, minorAllelePloidy);

        final double eventPenalty = EventPenalty.penalty(ploidyPenaltyFactor, majorAllelePloidy, minorAllelePloidy);
        final double deviationPenalty = (minorAllelePloidyDeviation + majorAllelePloidyDeviation) * observedBAF;

        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                .from(observedRegion)
                .fittedBAF(0)
                .fittedTumorCopyNumber(0)
                .tumorCopyNumber(impliedCopyNumber)
                .tumorBAF(impliedBAF)
                .refNormalisedCopyNumber(Doubles.replaceNaNWithZero(refNormalisedCopyNumber))
                .minorAlleleCopyNumberDeviation(minorAllelePloidyDeviation)
                .majorAlleleCopyNumberDeviation(majorAllelePloidyDeviation)
                .deviationPenalty(deviationPenalty)
                .eventPenalty(eventPenalty);

        return builder.build();
    }

    private double impliedBaf(final PurityAdjuster purityAdjuster, final String chromosome, final double copyNumber,
            final double observedBAF) {
        if (!cobaltChromosomes.contains(chromosome)) {
            return 1;
        }

        CobaltChromosome cobaltChromosome = cobaltChromosomes.get(chromosome);
        if (!cobaltChromosome.isNormal() || !cobaltChromosome.isDiploid()  || Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }

        return Doubles.lessOrEqual(observedBAF, ambiguousBaf)
                ? bafToMinimiseDeviation(purityAdjuster, chromosome, copyNumber)
                : purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF);
    }

    private double bafToMinimiseDeviation(final PurityAdjuster purityAdjuster, final String chromosome, double impliedCopyNumber) {
        final double minBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, 0.5)));
        final double maxBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, ambiguousBaf)));

        // Major Ploidy
        final double minBAFMajorAllelePloidy = minBAF * impliedCopyNumber;
        final double maxBAFMajorAllelePloidy = maxBAF * impliedCopyNumber;

        // Major Ploidy crosses whole number?
        final double minBAFMajorAllelePloidyCeil = Math.ceil(minBAFMajorAllelePloidy);
        if (!Doubles.equal(Math.signum(minBAFMajorAllelePloidyCeil - minBAFMajorAllelePloidy),
                Math.signum(minBAFMajorAllelePloidyCeil - maxBAFMajorAllelePloidy))) {
            return minBAFMajorAllelePloidyCeil / impliedCopyNumber;
        }

        // Minor Ploidy
        final double minBAFMinorAllelePloidy = impliedCopyNumber - minBAFMajorAllelePloidy;
        final double maxBAFMinorAllelePloidy = impliedCopyNumber - maxBAFMajorAllelePloidy;

        // Minor Ploidy crosses whole number?
        final double maxBAFMinorAllelePloidyCeil = Math.ceil(maxBAFMinorAllelePloidy);
        if (!Doubles.equal(Math.signum(maxBAFMinorAllelePloidyCeil - minBAFMinorAllelePloidy),
                Math.signum(maxBAFMinorAllelePloidyCeil - maxBAFMinorAllelePloidy))) {
            return 1 - maxBAFMinorAllelePloidyCeil / impliedCopyNumber;
        }

        double purity = purityAdjuster.purity();
        double normFactor = purityAdjuster.normFactor();

        // Minimise
        final double minBAFTotalDeviation =
                ploidyDeviation.majorAlleleDeviation(purity, normFactor, minBAFMajorAllelePloidy) + ploidyDeviation.minorAlleleDeviation(
                        purity,
                        normFactor,
                        minBAFMinorAllelePloidy);
        final double maxBAFTotalDeviation =
                ploidyDeviation.majorAlleleDeviation(purity, normFactor, maxBAFMajorAllelePloidy) + ploidyDeviation.minorAlleleDeviation(
                        purity,
                        normFactor,
                        maxBAFMinorAllelePloidy);
        return Doubles.lessThan(minBAFTotalDeviation, maxBAFTotalDeviation) ? 0.5 : ambiguousBaf;
    }
}