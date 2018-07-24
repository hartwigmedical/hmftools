package com.hartwig.hmftools.common.purple.region;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public class FittedRegionFactoryV2 {

    private static final double AMBIGUOUS_BAF = 0.542;

    private final Gender gender;

    public FittedRegionFactoryV2(final Gender gender) {
        this.gender = gender;
    }

    @NotNull
    public List<FittedRegion> fitRegion(final double purity, final double normFactor,
            @NotNull final Collection<ObservedRegion> observedRegions) {

        final Predicate<ObservedRegion> valid = observedRegion -> gender == Gender.MALE || !observedRegion.chromosome().equals("Y");
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @NotNull
    public FittedRegion fitRegion(final double purity, final double normFactor, final @NotNull ObservedRegion observedRegion) {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(gender, purity, normFactor);

        double observedTumorRatio = observedRegion.observedTumorRatio();
        double impliedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double observedBAF = observedRegion.observedBAF();
        double impliedBAF = impliedBaf(purityAdjuster, observedRegion.chromosome(), impliedCopyNumber, observedBAF);

        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());

        double majorAllelePloidy = impliedBAF * impliedCopyNumber;
        double majorAllelePloidyDeviation = majorDeviation(majorAllelePloidy);
        double minorAllelePloidy = impliedCopyNumber - majorAllelePloidy;
        double minorAllelePloidyDeviation = minorDeviation(minorAllelePloidy);

        double ploidyPenalty = PloidyPenalty.penalty((int) Math.round(impliedCopyNumber), (int) Math.round(majorAllelePloidy));
        double totalDeviation = minorAllelePloidyDeviation + majorAllelePloidyDeviation;

        ImmutableFittedRegion.Builder builder = ImmutableFittedRegion.builder()
                .from(observedRegion)
                .segmentBAF(0)
                .segmentTumorCopyNumber(0)
                .tumorCopyNumber(impliedCopyNumber)
                .tumorBAF(impliedBAF)
                .refNormalisedCopyNumber(Doubles.replaceNaNWithZero(refNormalisedCopyNumber))
                .modelBAF(0)
                .modelPloidy((int)Math.round(majorAllelePloidy))
                .modelTumorRatio(0)
                .bafDeviation(majorAllelePloidyDeviation)
                .cnvDeviation(minorAllelePloidyDeviation)
                .deviation(totalDeviation)
                .ploidyPenalty(ploidyPenalty);

        return builder.build();
    }

    public double impliedBaf(final PurityAdjuster purityAdjuster, final String chromosome, final double copyNumber,
            final double observedBAF) {
        boolean isHomologous = HumanChromosome.fromString(chromosome).isDiploid(gender);

        if (!isHomologous || Doubles.lessOrEqual(copyNumber, 1)) {
            return 1;
        }

        return Doubles.lessOrEqual(observedBAF, AMBIGUOUS_BAF)
                ? bafToMinimiseDeviation(purityAdjuster, chromosome, copyNumber)
                : purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF);

    }

    @VisibleForTesting
    static double bafToMinimiseDeviation(final PurityAdjuster purityAdjuster, final String chromosome, double impliedCopyNumber) {

        final double minBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, 0.5)));
        final double maxBAF =
                Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, AMBIGUOUS_BAF)));

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

        // Minimise
        final double minBAFTotalDeviation = majorDeviation(minBAFMajorAllelePloidy) + minorDeviation(minBAFMinorAllelePloidy);
        final double maxBAFTotalDeviation = majorDeviation(maxBAFMajorAllelePloidy) + minorDeviation(maxBAFMinorAllelePloidy);
        return Doubles.lessThan(minBAFTotalDeviation, maxBAFTotalDeviation) ? 0.5 : AMBIGUOUS_BAF;
    }

    private static double majorDeviation(double majorAllelePloidy) {
        return deviationFromWholeNumber(1, majorAllelePloidy);
    }

    private static double minorDeviation(double minorAllelePloidy) {
        return deviationFromWholeNumber(0, minorAllelePloidy);
    }

    private static double deviationFromWholeNumber(int min, double impliedPloidy) {
        long wholeNumber = Math.max(min, Math.round(impliedPloidy));
        return Math.abs(impliedPloidy - wholeNumber);
    }

}