package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.HIGH_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_2_PLUS_FRAGS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_2_PLUS_FRAG_PERC;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_FRAG_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityData.NO_RESULT;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.cappedPurity;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.estimatedPurity;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class LowCountModel extends ClonalityModel
{
    public LowCountModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        super(config, resultsWriter, sample,  variants);
    }

    public static List<SomaticVariant> filterVariants(
            final String sampleId, final FragmentTotals fragmentTotals, final List<SomaticVariant> variants, double medianVcn)
    {
        // We should exclude variants from LOW_COUNT which are not close to the normal copy number profile or depth of the sample.
        // ie. if VCN > 2x median VCN or if sampleDP > 2x wAD
        double vcnThreshold = 2 * medianVcn;
        double sampleDpThreshold = 2 * fragmentTotals.weightedSampleDepth();

        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        for(SomaticVariant variant : variants)
        {
            if(variant.VariantCopyNumber > vcnThreshold)
                continue;

            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(sampleFragData.Depth > sampleDpThreshold)
                continue;

            filteredVariants.add(variant);
        }

        return filteredVariants;
    }

    public static boolean canUseModel(final String sampleId, final FragmentTotals fragmentTotals, final List<SomaticVariant> variants)
    {
        int twoPlusFragVariantCount = 0;

        for(SomaticVariant variant : variants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(sampleFragData.AlleleCount >= 2)
                ++twoPlusFragVariantCount;
        }

        if(twoPlusFragVariantCount < LOW_COUNT_MODEL_MIN_FRAG_VARIANTS)
            return false;

        if(fragmentTotals.weightedSampleDepth() >= LOW_COUNT_MODEL_MIN_AVG_DEPTH)
            return false;

        double twoPlusPercent = twoPlusFragVariantCount / (double)variants.size();

        return twoPlusFragVariantCount >= LOW_COUNT_MODEL_MIN_2_PLUS_FRAGS && twoPlusPercent >= LOW_COUNT_MODEL_MIN_2_PLUS_FRAG_PERC;
    }

    @Override
    public void calculate(final String sampleId, final FragmentTotals fragmentTotals, final PurityCalcData purityCalcData, double noiseRate)
    {
        ClonalityResult clonalityResult = calcClonalityData(sampleId, fragmentTotals);

        if(clonalityResult == null)
            return;

        String clonalityInfo = format("Dropout(rate=%.2f weightedCN=%.2f weightedVCN=%.2f)",
                clonalityResult.dropoutRate(), clonalityResult.weightedCn(), clonalityResult.weightedVcn());

        purityCalcData.Clonality = new ClonalityData(
                ClonalityMethod.LOW_COUNT, clonalityResult.dropoutVaf(), clonalityResult.dropoutVafLow(), clonalityResult.dropoutVafHigh(),
                clonalityResult.variantCount(), clonalityInfo);

        double lowRatio = purityCalcData.PurityEstimate > 0 ? purityCalcData.PurityRangeLow / purityCalcData.PurityEstimate : 1;
        double highRatio = purityCalcData.PurityEstimate > 0 ? purityCalcData.PurityRangeHigh / purityCalcData.PurityEstimate : 1;
        double INVALID_PURITY = -1;

        double purityEstimate = cappedPurity(estimatedPurity(
                purityCalcData.Clonality.Vaf, noiseRate, clonalityResult.weightedVcn, clonalityResult.weightedCn), INVALID_PURITY);

        if(purityEstimate != INVALID_PURITY)
        {
            purityCalcData.PurityEstimate = purityEstimate;

            purityCalcData.PurityRangeLow = cappedPurity(estimatedPurity(
                    purityCalcData.Clonality.VafLow, noiseRate, fragmentTotals), purityEstimate);

            purityCalcData.PurityRangeHigh = cappedPurity(estimatedPurity(
                    purityCalcData.Clonality.VafHigh, noiseRate, fragmentTotals), purityEstimate);
        }

        purityCalcData.PurityRangeLow = min(purityCalcData.PurityRangeLow * lowRatio, 1);
        purityCalcData.PurityRangeHigh = min(purityCalcData.PurityRangeHigh * highRatio, 1);
    }

    private record ClonalityResult(
            double weightedVcn, double weightedCn, double dropoutRate, double dropoutVaf,
            double dropoutVafLow, double dropoutVafHigh, int variantCount) {}

    private ClonalityResult calcClonalityData(final String sampleId, final FragmentTotals fragmentTotals)
    {
        if(fragmentTotals.oneFragmentCount() == 0)
            return null;

        double estimateVaf = fragmentTotals.adjSampleVaf();

        if(estimateVaf == 0)
            return null;

        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        int observedFrag1 = 0;
        int observedFrag2Plus = 0;

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(useVariant(variant, sampleFragData) && sampleFragData.Depth > 1)
            {
                filteredVariants.add(variant);

                if(sampleFragData.AlleleCount == 1)
                    ++observedFrag1;
                else if(sampleFragData.AlleleCount >= 2)
                    ++observedFrag2Plus;
            }
        }

        List<SimulatedVafCalcs> simulatedVafCalcs = Lists.newArrayList();

        // now test each simulated dropout rate and VAF
        for(double dropoutRate = 0; dropoutRate < 0.95; dropoutRate += PurityConstants.DROPOUT_RATE_INCREMENT)
        {
            double simulatedVaf = estimateVaf / (1 - dropoutRate);

            if(simulatedVaf >= 1)
                break;

            double probTotalFrag0 = 0;
            double probTotalFrag1 = 0;
            double probTotalFrag2Plus = 0;
            double depthToCopyNumberFactor;

            double depthToCopyNumberFactorSumNon0 = 0;
            double weightedVcnSumNon0 = 0;
            double weightedCnSumNon0 = 0;
            double depthToCopyNumberFactorSum0 = 0;
            double weightedVcnSum0 = 0;
            double weightedCnSum0 = 0;

            for(SomaticVariant variant : filteredVariants)
            {
                GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

                BinomialDistribution binomial = new BinomialDistribution(sampleFragData.Depth, simulatedVaf);

                double probFrag0 = binomial.probability(0);
                double probFrag1 = binomial.probability(1);
                double probFrag2Plus = 1 - binomial.cumulativeProbability(1);

                probTotalFrag0 += probFrag0;
                probTotalFrag1 += probFrag1;
                probTotalFrag2Plus += probFrag2Plus;
                depthToCopyNumberFactor = sampleFragData.Depth / max(variant.CopyNumber, 1); // should link this to the definition in FragmentTotals somehow

                if(sampleFragData.AlleleCount > 0) // non-dropout variants
                {
                    depthToCopyNumberFactorSumNon0 += depthToCopyNumberFactor;
                    weightedVcnSumNon0 += variant.VariantCopyNumber * depthToCopyNumberFactor;
                    weightedCnSumNon0 += variant.CopyNumber * depthToCopyNumberFactor;
                }
                else
                {
                    depthToCopyNumberFactorSum0 += depthToCopyNumberFactor;
                    weightedVcnSum0 += variant.VariantCopyNumber * depthToCopyNumberFactor;
                    weightedCnSum0 += variant.CopyNumber * depthToCopyNumberFactor;
                }
            }

            double zeroWeighting = probTotalFrag0 / (probTotalFrag0 + probTotalFrag1 + probTotalFrag2Plus);

            double weightedVcn = zeroWeighting * (weightedVcnSum0 / depthToCopyNumberFactorSum0)
                    + (1 - zeroWeighting) * (weightedVcnSumNon0 / depthToCopyNumberFactorSumNon0);

            double weightedCn = zeroWeighting * (weightedCnSum0 / depthToCopyNumberFactorSum0)
                    + (1 - zeroWeighting) * (weightedCnSumNon0 / depthToCopyNumberFactorSumNon0);

            SimulatedVafCalcs simVafCalcs = new SimulatedVafCalcs(
                    dropoutRate, simulatedVaf, probTotalFrag1, probTotalFrag2Plus, weightedVcn, weightedCn);
            simulatedVafCalcs.add(simVafCalcs);
        }

        // now find the closest ratio to the observed ratio
        double observedRatio = observedFrag2Plus / (double)observedFrag1;
        DropoutVaf dropoutVaf = findVafRatio(observedRatio, simulatedVafCalcs);

        if(dropoutVaf.Dropout == 0)
            return null;

        // LOW_COUNT_PROBABILITY
        double lowFragmentCount = max(calcPoissonNoiseValue(observedFrag2Plus, HIGH_PROBABILITY), 2);
        double highFragmentCount = calcPoissonNoiseValue(observedFrag2Plus, LOW_PROBABILITY);

        double lowObservedRatio = lowFragmentCount / observedFrag1;
        DropoutVaf dropoutVafLow = findVafRatio(lowObservedRatio, simulatedVafCalcs);

        double highObservedRatio = highFragmentCount / observedFrag1;
        DropoutVaf dropoutVafHigh = findVafRatio(highObservedRatio, simulatedVafCalcs);

        CT_LOGGER.debug(format("sample(%s) low-count model: obsRatio(%.2f 1=%d 2+=%d) vaf((%.6f low=%.6f high=%.6f)",
                sampleId, observedRatio, observedFrag1, observedFrag2Plus, dropoutVaf.VAF, dropoutVafLow.VAF, dropoutVafHigh.VAF));

        /*
         return new ClonalityData(
        ClonalityMethod.LOW_COUNT, dropoutVaf.VAF, dropoutVafLow.VAF, dropoutVafHigh.VAF,
        observedFrag1 + observedFrag2Plus, dropoutVaf.Dropout,0, 0, 0,
        dropoutVaf.Vcn, dropoutVaf.Cn);
         */

        return new ClonalityResult(
                dropoutVaf.Vcn, dropoutVaf.Cn, dropoutVaf.Dropout, dropoutVaf.VAF, dropoutVafLow.VAF, dropoutVafHigh.VAF,
                observedFrag1 + observedFrag2Plus);
    }

    private DropoutVaf findVafRatio(double observedRatio, final List<SimulatedVafCalcs> simulatedVafCalcs)
    {
        SimulatedVafCalcs closestRatioUp = null;
        SimulatedVafCalcs closestRatioDown = null;

        for(SimulatedVafCalcs simVafCalcs : simulatedVafCalcs)
        {
            double probRatio = simVafCalcs.fragmentRatio();

            if(probRatio > observedRatio)
            {
                if(closestRatioUp == null || probRatio < closestRatioUp.fragmentRatio())
                    closestRatioUp = simVafCalcs;
            }
            else
            {
                if(closestRatioDown == null || probRatio > closestRatioDown.fragmentRatio())
                    closestRatioDown = simVafCalcs;
            }
        }

        double calcVaf = 0;
        double calcDropout = 0;
        double calcVcn = 0;
        double calcCn = 0;

        if(closestRatioUp != null && closestRatioDown != null)
        {
            double upperFraction = (observedRatio - closestRatioDown.fragmentRatio())
                    / (closestRatioUp.fragmentRatio() - closestRatioDown.fragmentRatio());

            calcVaf = upperFraction * closestRatioUp.SimulatedVaf + (1 - upperFraction) * closestRatioDown.SimulatedVaf;
            calcDropout = upperFraction * closestRatioUp.DropoutRate + (1 - upperFraction) * closestRatioDown.DropoutRate;
            calcVcn = upperFraction * closestRatioUp.WeightedVcn + (1 - upperFraction) * closestRatioDown.WeightedVcn;
            calcCn = upperFraction * closestRatioUp.WeightedCn + (1 - upperFraction) * closestRatioDown.WeightedCn;

        }
        else if(closestRatioUp != null)
        {
            calcVaf = closestRatioUp.SimulatedVaf;
            calcDropout = closestRatioUp.DropoutRate;
            calcVcn = closestRatioUp.WeightedVcn;
            calcCn = closestRatioUp.WeightedCn;
        }
        else
        {
            calcVaf = closestRatioDown.SimulatedVaf;
            calcDropout = closestRatioDown.DropoutRate;
            calcVcn = closestRatioDown.WeightedVcn;
            calcCn = closestRatioDown.WeightedCn;
        }

        return new DropoutVaf(calcVaf, calcDropout, calcVcn, calcCn);
    }

    private class DropoutVaf
    {
        public final double VAF;
        public final double Dropout;
        public final double Vcn;
        public final double Cn;

        public DropoutVaf(final double vaf, final double dropout, double vcn, double cn)
        {
            VAF = vaf;
            Dropout = dropout;
            Vcn = vcn;
            Cn = cn;
        }
    }

    private class SimulatedVafCalcs
    {
        public final double DropoutRate;
        public final double SimulatedVaf;
        public final double ProbabilityTotalFrag1;
        public final double ProbabilityTotalFrag2Plus;
        public final double WeightedVcn;
        public final double WeightedCn;

        public SimulatedVafCalcs(
                final double dropoutRate, final double simulatedVaf, final double probabilityTotalFrag1,
                final double probabilityTotalFrag2Plus, double weightedVcn, double weightedCn)
        {
            DropoutRate = dropoutRate;
            SimulatedVaf = simulatedVaf;
            ProbabilityTotalFrag1 = probabilityTotalFrag1;
            ProbabilityTotalFrag2Plus = probabilityTotalFrag2Plus;
            WeightedVcn = weightedVcn;
            WeightedCn = weightedCn;
        }

        public double fragmentRatio() { return ProbabilityTotalFrag1 > 0 ? ProbabilityTotalFrag2Plus / ProbabilityTotalFrag1 : 0; }

        public String toString() { return format("simVaf(%.4f) doRate(%.2f) ratio(%.3f)", SimulatedVaf, DropoutRate, fragmentRatio()); }
    }
}
