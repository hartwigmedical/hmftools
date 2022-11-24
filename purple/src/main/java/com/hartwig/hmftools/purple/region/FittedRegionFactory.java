package com.hartwig.hmftools.purple.region;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class FittedRegionFactory
{
    private final double mAmbiguousBaf;
    private final double mPloidyPenaltyFactor;
    private final PloidyDeviation mPloidyDeviation;
    private final CobaltChromosomes mCobaltChromosomes;

    public FittedRegionFactory(final CobaltChromosomes cobaltChromosomes, final int averageReadDepth, double ploidyPenaltyFactor,
            double ploidyPenaltyStandardDeviation, double ploidyPenaltyMinStandardDeviationPerPloidy,
            final double majorAlleleSubOnePenaltyMultiplier, final double majorAlleleSubOneAdditionalPenalty,
            final double baselineDeviation)
    {
        mCobaltChromosomes = cobaltChromosomes;
        mPloidyPenaltyFactor = ploidyPenaltyFactor;
        mPloidyDeviation = new PloidyDeviation(ploidyPenaltyStandardDeviation,
                ploidyPenaltyMinStandardDeviationPerPloidy,
                majorAlleleSubOnePenaltyMultiplier,
                majorAlleleSubOneAdditionalPenalty,
                baselineDeviation);
        mAmbiguousBaf = ExpectedBAF.expectedBAF(averageReadDepth);
    }

    public List<ObservedRegion> fitRegion(double purity, double normFactor, @NotNull final Collection<ObservedRegion> observedRegions)
    {
        final Predicate<ObservedRegion> valid = observedRegion -> isAllowedRegion(mCobaltChromosomes, observedRegion);
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean isAllowedRegion(@NotNull final CobaltChromosomes cobaltChromosomes, final GenomeRegion region)
    {
        return cobaltChromosomes.contains(region.chromosome());
    }

    public ObservedRegion fitRegion(final double purity, final double normFactor, final ObservedRegion observedRegion)
    {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purity, normFactor, mCobaltChromosomes.chromosomes(), mCobaltChromosomes.gender());

        double observedTumorRatio = observedRegion.observedTumorRatio();
        double impliedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double observedBAF = observedRegion.observedBAF();
        double impliedBAF = impliedBaf(purityAdjuster, observedRegion.chromosome(), impliedCopyNumber, observedBAF);

        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());

        double majorAllelePloidy = impliedBAF * impliedCopyNumber;
        double minorAllelePloidy = impliedCopyNumber - majorAllelePloidy;

        double majorAllelePloidyDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAllelePloidy);
        double minorAllelePloidyDeviation = mPloidyDeviation.minorAlleleDeviation(purity, normFactor, minorAllelePloidy);

        final double eventPenalty = EventPenalty.penalty(mPloidyPenaltyFactor, majorAllelePloidy, minorAllelePloidy);
        final double deviationPenalty = (minorAllelePloidyDeviation + majorAllelePloidyDeviation) * observedBAF;

        ObservedRegion fittedRegion = ObservedRegion.from(observedRegion);

        fittedRegion.setTumorCopyNumber(impliedCopyNumber);

        fittedRegion.setTumorBAF(impliedBAF);
        fittedRegion.setRefNormalisedCopyNumber(Doubles.replaceNaNWithZero(refNormalisedCopyNumber));
        fittedRegion.setMinorAlleleCopyNumberDeviation(minorAllelePloidyDeviation);
        fittedRegion.setMajorAlleleCopyNumberDeviation(majorAllelePloidyDeviation);
        fittedRegion.setDeviationPenalty(deviationPenalty);
        fittedRegion.setEventPenalty(eventPenalty);

        return fittedRegion;
    }

    private static final double MIN_CN_THRESHOLD = 0.1;
    // private static final double MIN_CN_THRESHOLD = 1;

    private double impliedBaf(final PurityAdjuster purityAdjuster, final String chromosome, final double copyNumber,
            final double observedBAF)
    {
        if(!mCobaltChromosomes.contains(chromosome))
        {
            return 1;
        }

        CobaltChromosome cobaltChromosome = mCobaltChromosomes.get(chromosome);
        if(!cobaltChromosome.isNormal() || !cobaltChromosome.isDiploid() || Doubles.lessOrEqual(copyNumber, MIN_CN_THRESHOLD))
        {
            return 1;
        }

        return Doubles.lessOrEqual(observedBAF, mAmbiguousBaf)
                ? bafToMinimiseDeviation(purityAdjuster, chromosome, copyNumber)
                : purityAdjuster.purityAdjustedBAFSimple(chromosome, copyNumber, observedBAF);
    }

    private double bafToMinimiseDeviation(final PurityAdjuster purityAdjuster, final String chromosome, double impliedCopyNumber)
    {
        final double minBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, 0.5)));
        final double maxBAF = Math.max(0, Math.min(1, purityAdjuster.purityAdjustedBAFSimple(chromosome, impliedCopyNumber, mAmbiguousBaf)));

        // Major Ploidy
        final double minBAFMajorAllelePloidy = minBAF * impliedCopyNumber;
        final double maxBAFMajorAllelePloidy = maxBAF * impliedCopyNumber;

        // Major Ploidy crosses whole number?
        final double minBAFMajorAllelePloidyCeil = Math.ceil(minBAFMajorAllelePloidy);
        if(!Doubles.equal(Math.signum(minBAFMajorAllelePloidyCeil - minBAFMajorAllelePloidy),
                Math.signum(minBAFMajorAllelePloidyCeil - maxBAFMajorAllelePloidy)))
        {
            return minBAFMajorAllelePloidyCeil / impliedCopyNumber;
        }

        // Minor Ploidy
        final double minBAFMinorAllelePloidy = impliedCopyNumber - minBAFMajorAllelePloidy;
        final double maxBAFMinorAllelePloidy = impliedCopyNumber - maxBAFMajorAllelePloidy;

        // Minor Ploidy crosses whole number?
        final double maxBAFMinorAllelePloidyCeil = Math.ceil(maxBAFMinorAllelePloidy);
        if(!Doubles.equal(Math.signum(maxBAFMinorAllelePloidyCeil - minBAFMinorAllelePloidy),
                Math.signum(maxBAFMinorAllelePloidyCeil - maxBAFMinorAllelePloidy)))
        {
            return 1 - maxBAFMinorAllelePloidyCeil / impliedCopyNumber;
        }

        double purity = purityAdjuster.purity();
        double normFactor = purityAdjuster.normFactor();

        // Minimise
        final double minBAFTotalDeviation =
                mPloidyDeviation.majorAlleleDeviation(purity, normFactor, minBAFMajorAllelePloidy) + mPloidyDeviation.minorAlleleDeviation(
                        purity,
                        normFactor,
                        minBAFMinorAllelePloidy);
        final double maxBAFTotalDeviation =
                mPloidyDeviation.majorAlleleDeviation(purity, normFactor, maxBAFMajorAllelePloidy) + mPloidyDeviation.minorAlleleDeviation(
                        purity,
                        normFactor,
                        maxBAFMinorAllelePloidy);
        return Doubles.lessThan(minBAFTotalDeviation, maxBAFTotalDeviation) ? 0.5 : mAmbiguousBaf;
    }
}