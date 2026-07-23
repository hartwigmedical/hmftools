package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
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

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConstants;

import org.apache.commons.math3.distribution.PoissonDistribution;

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
    public ClonalityData calculate(final String sampleId, final FragmentTotals fragmentTotals, final PurityCalcData purityCalcData)
    {
        double estimateVaf = fragmentTotals.adjSampleVaf();

        if(estimateVaf == 0)
            return NO_RESULT;

        List<SomaticVariant> filteredVariants = Lists.newArrayList();

        int observedFrag1 = fragmentTotals.sampleOneFragmentCount();
        int observedFrag2Plus = fragmentTotals.sampleTwoPlusCount();

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(useVariant(variant, sampleFragData) && sampleFragData.Depth > 0)
            {
                filteredVariants.add(variant);

                if(sampleFragData.AlleleCount == 1)
                    ++observedFrag1;
                else if(sampleFragData.AlleleCount >= 2)
                    ++observedFrag2Plus;
            }
        }

        if(fragmentTotals.sampleOneFragmentCount() == 0)
            return NO_RESULT;

        List<SimulatedVafCalcs> simulatedVafCalcs = Lists.newArrayList();

        // now test each simulated dropout rate and VAF
        for(double dropoutRate = 0; dropoutRate < 0.95; dropoutRate += PurityConstants.DROPOUT_RATE_INCREMENT)
        {
            double simulatedVaf = estimateVaf / (1 - dropoutRate);

            if(simulatedVaf >= 1)
                break;

            double probTotalFrag1 = 0;
            double probTotalFrag2Plus = 0;

            for(SomaticVariant variant : filteredVariants)
            {
                GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

                double expectedAD = sampleFragData.Depth * simulatedVaf;

                PoissonDistribution poisson = new PoissonDistribution(expectedAD);

                double probFrag1 = poisson.probability(1);
                double probFrag2Plus = 1 - poisson.cumulativeProbability(1);

                probTotalFrag1 += probFrag1;
                probTotalFrag2Plus += probFrag2Plus;
            }

            SimulatedVafCalcs simVafCalcs = new SimulatedVafCalcs(dropoutRate, simulatedVaf, probTotalFrag1, probTotalFrag2Plus);
            simulatedVafCalcs.add(simVafCalcs);
        }

        // now find the closest ratio to the observed ratio
        double observedRatio = observedFrag2Plus / (double)observedFrag1;
        DropoutVaf dropoutVaf = findVafRatio(observedRatio, simulatedVafCalcs);

        if(dropoutVaf.Dropout == 0)
            return NO_RESULT;

        // LOW_COUNT_PROBABILITY
        double lowFragmentCount = max(calcPoissonNoiseValue(observedFrag2Plus, HIGH_PROBABILITY), 2);
        double highFragmentCount = calcPoissonNoiseValue(observedFrag2Plus, LOW_PROBABILITY);

        double lowObservedRatio = lowFragmentCount / observedFrag1;
        DropoutVaf dropoutVafLow = findVafRatio(lowObservedRatio, simulatedVafCalcs);

        double highObservedRatio = highFragmentCount / observedFrag1;
        DropoutVaf dropoutVafHigh = findVafRatio(highObservedRatio, simulatedVafCalcs);

        CT_LOGGER.debug(format("sample(%s) low-count model: obsRatio(%.2f 1=%d 2+=%d) vaf((%.6f low=%.6f high=%.6f)",
                sampleId, observedRatio, observedFrag1, observedFrag2Plus, dropoutVaf.VAF, dropoutVafLow.VAF, dropoutVafHigh.VAF));

        return new ClonalityData(
                ClonalityMethod.LOW_COUNT, dropoutVaf.VAF, dropoutVafLow.VAF, dropoutVafHigh.VAF,
                observedFrag1 + observedFrag2Plus, dropoutVaf.Dropout,0, 0, 0);
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

        if(closestRatioUp != null && closestRatioDown != null)
        {
            double upperFraction = (observedRatio - closestRatioDown.fragmentRatio())
                    / (closestRatioUp.fragmentRatio() - closestRatioDown.fragmentRatio());

            calcVaf = upperFraction * closestRatioUp.SimulatedVaf + (1 - upperFraction) * closestRatioDown.SimulatedVaf;
            calcDropout = upperFraction * closestRatioUp.DropoutRate + (1 - upperFraction) * closestRatioDown.DropoutRate;
        }
        else if(closestRatioUp != null)
        {
            calcVaf = closestRatioUp.SimulatedVaf;
            calcDropout = closestRatioUp.DropoutRate;
        }
        else
        {
            calcVaf = closestRatioDown.SimulatedVaf;
            calcDropout = closestRatioDown.DropoutRate;
        }

        return new DropoutVaf(calcVaf, calcDropout);
    }

    private class DropoutVaf
    {
        public final double VAF;
        public final double Dropout;

        public DropoutVaf(final double vaf, final double dropout)
        {
            VAF = vaf;
            Dropout = dropout;
        }
    }

    private class SimulatedVafCalcs
    {
        public final double DropoutRate;
        public final double SimulatedVaf;
        public final double ProbabilityTotalFrag1;
        public final double ProbabilityTotalFrag2Plus;

        public SimulatedVafCalcs(
                final double dropoutRate, final double simulatedVaf, final double probabilityTotalFrag1,
                final double probabilityTotalFrag2Plus)
        {
            DropoutRate = dropoutRate;
            SimulatedVaf = simulatedVaf;
            ProbabilityTotalFrag1 = probabilityTotalFrag1;
            ProbabilityTotalFrag2Plus = probabilityTotalFrag2Plus;
        }

        public double fragmentRatio() { return ProbabilityTotalFrag1 > 0 ? ProbabilityTotalFrag2Plus / ProbabilityTotalFrag1 : 0; }

        public String toString() { return format("simVaf(%.4f) doRate(%.2f) ratio(%.3f)", SimulatedVaf, DropoutRate, fragmentRatio()); }
    }
}
