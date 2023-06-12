package com.hartwig.hmftools.purple.purity;

import java.util.Collection;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.region.PloidyDeviation;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
import com.hartwig.hmftools.common.utils.Doubles;

public class RegionFitCalculator
{
    private final double mAmbiguousBaf;
    private final double mPloidyPenaltyFactor;
    private final PloidyDeviation mPloidyDeviation;
    private final CobaltChromosomes mCobaltChromosomes;

    public RegionFitCalculator(final CobaltChromosomes cobaltChromosomes, final int averageReadDepth, double ploidyPenaltyFactor,
            double ploidyPenaltyStandardDeviation, double ploidyPenaltyMinStandardDeviationPerPloidy,
            final double majorAlleleSubOnePenaltyMultiplier, final double majorAlleleSubOneAdditionalPenalty,
            final double baselineDeviation)
    {
        mCobaltChromosomes = cobaltChromosomes;
        mPloidyPenaltyFactor = ploidyPenaltyFactor;

        mPloidyDeviation = new PloidyDeviation(
                ploidyPenaltyStandardDeviation, ploidyPenaltyMinStandardDeviationPerPloidy,
                majorAlleleSubOnePenaltyMultiplier, majorAlleleSubOneAdditionalPenalty,
                baselineDeviation);

        // mPloidyDeviation.setUseCache();

        mAmbiguousBaf = ExpectedBAF.expectedBAF(averageReadDepth);
    }

    public List<ObservedRegion> fitRegion(double purity, double normFactor, final Collection<ObservedRegion> observedRegions)
    {
        final Predicate<ObservedRegion> valid = observedRegion -> isAllowedRegion(mCobaltChromosomes, observedRegion);
        return observedRegions.stream().filter(valid).map(x -> fitRegion(purity, normFactor, x)).collect(Collectors.toList());
    }

    @VisibleForTesting
    public static boolean isAllowedRegion(final CobaltChromosomes cobaltChromosomes, final GenomeRegion region)
    {
        return cobaltChromosomes.contains(region.chromosome());
    }

    public ObservedRegion fitRegion(final double purity, final double normFactor, final ObservedRegion observedRegion)
    {
        RegionFitCalcs regionFitCalcs = calculateRegionFit(purity, normFactor, observedRegion);

        ObservedRegion fittedRegion = ObservedRegion.from(observedRegion);

        fittedRegion.setTumorCopyNumber(regionFitCalcs.TumorCopyNumber);
        fittedRegion.setTumorBAF(regionFitCalcs.TumorBAF);
        fittedRegion.setRefNormalisedCopyNumber(regionFitCalcs.RefNormalisedCopyNumber);
        fittedRegion.setMinorAlleleCopyNumberDeviation(regionFitCalcs.MinorAlleleCopyNumberDeviation);
        fittedRegion.setMajorAlleleCopyNumberDeviation(regionFitCalcs.MajorAlleleCopyNumberDeviation);
        fittedRegion.setDeviationPenalty(regionFitCalcs.DeviationPenalty);
        fittedRegion.setEventPenalty(regionFitCalcs.EventPenalty);

        return fittedRegion;
    }

    public RegionFitCalcs calculateRegionFit(final double purity, final double normFactor, final ObservedRegion observedRegion)
    {
        final PurityAdjuster purityAdjuster = new PurityAdjuster(purity, normFactor, mCobaltChromosomes);

        double observedTumorRatio = observedRegion.observedTumorRatio();
        double impliedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedRegion.chromosome(), observedTumorRatio);
        double observedBAF = observedRegion.observedBAF();
        double impliedBAF = impliedBaf(purityAdjuster, observedRegion.chromosome(), impliedCopyNumber, observedBAF);

        double refNormalisedCopyNumber = purityAdjuster.purityAdjustedCopyNumber(observedTumorRatio, observedRegion.observedNormalRatio());

        double majorAllelePloidy = impliedBAF * impliedCopyNumber;
        double minorAllelePloidy = impliedCopyNumber - majorAllelePloidy;

        double majorAllelePloidyDeviation = mPloidyDeviation.majorAlleleDeviation(purity, normFactor, majorAllelePloidy);
        double minorAllelePloidyDeviation = mPloidyDeviation.minorAlleleDeviation(purity, normFactor, minorAllelePloidy);

        double eventPenalty = calculateEventPenalty(mPloidyPenaltyFactor, majorAllelePloidy, minorAllelePloidy);
        double deviationPenalty = (minorAllelePloidyDeviation + majorAllelePloidyDeviation) * observedBAF;

        return new RegionFitCalcs(
                impliedCopyNumber, impliedBAF, Doubles.replaceNaNWithZero(refNormalisedCopyNumber),
                minorAllelePloidyDeviation, majorAllelePloidyDeviation, eventPenalty, deviationPenalty);
    }

    public static double calculateEventPenalty(double eventPenaltyFactor, double majorAllele, double minorAllele)
    {
        double wholeGenomeDoublingDistance = 1 + (Math.abs(majorAllele - 2)) + (Math.abs(minorAllele - 2));
        double singleEventDistance = (Math.abs(majorAllele - 1)) + (Math.abs(minorAllele - 1));

        return 1 + eventPenaltyFactor * Math.min(singleEventDistance, wholeGenomeDoublingDistance);
    }

    private static final double MIN_CN_THRESHOLD = 0.1;

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