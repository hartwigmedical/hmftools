package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.HIGH_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SNV_QUAL_THRESHOLDS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SYNTHETIC_TUMOR_VAF;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;
import static com.hartwig.hmftools.wisp.purity.variant.BqrAdjustment.hasVariantContext;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.convertVafToPurity;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.isRecomputed;
import static com.hartwig.hmftools.wisp.purity.variant.LowCountModel.filterVariants;
import static com.hartwig.hmftools.wisp.purity.variant.PurityCalcData.CALC_NO_SET;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.calcLimitOfDetection;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.estimatedProbability;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.estimatedPurity;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult.INVALID_RESULT;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class SomaticPurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final SampleData mSample;
    private final BqrAdjustment mBqrAdjustment;

    public SomaticPurityEstimator(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
        mBqrAdjustment = new BqrAdjustment(mConfig);
    }

    public SomaticPurityResult calculatePurity(
            final String sampleId, final PurityContext purityContext, final List<SomaticVariant> variants,
            final int totalVariantCount, final int chipVariants)
    {
        FragmentTotals fragmentTotals = new FragmentTotals();

        int sampleDualDP = 0;
        int sampleDualAD = 0;

        UmiTypeCounts umiTypeCounts = new UmiTypeCounts();

        for(SomaticVariant variant : variants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            fragmentTotals.addVariantData(
                    variant.CopyNumber, variant.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, sampleFragData.QualTotal);

            umiTypeCounts.add(sampleFragData.UmiCounts);
            sampleDualDP += sampleFragData.UmiCounts.TotalDual;
            sampleDualAD += sampleFragData.UmiCounts.AlleleDual;
        }

        if(fragmentTotals.sampleDepthTotal() == 0)
            return INVALID_RESULT;

        if(!mConfig.SkipBqr)
            mBqrAdjustment.loadBqrData(sampleId);

        if(mConfig.hasSyntheticTumor())
        {
            fragmentTotals.setTumorVafOverride(SYNTHETIC_TUMOR_VAF);
        }

        PurityCalcData purityCalcData = new PurityCalcData();

        // firstly estimate raw purity without consideration of clonal peaks
        double noiseRate = mConfig.noiseRate(false);

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 99% confidence of a tumor presence
        if(!mConfig.SkipBqr)
        {
            if(!mBqrAdjustment.hasValidData())
                return INVALID_RESULT;

            for(int threshold : SNV_QUAL_THRESHOLDS)
            {
                FragmentTotals bqrThresholdFragTotals = calculateThresholdValues(sampleId, purityCalcData, variants, threshold);

                if(bqrThresholdFragTotals != null)
                    fragmentTotals = bqrThresholdFragTotals;
            }

            noiseRate = purityCalcData.ErrorRate;
        }
        else
        {
            purityCalcData.RawPurityEstimate = estimatedPurity(fragmentTotals.rawSampleVaf(), noiseRate, fragmentTotals);
            purityCalcData.Probability = estimatedProbability(fragmentTotals, noiseRate);
            purityCalcData.LodPurityEstimate = calcLimitOfDetection(fragmentTotals, noiseRate);
            purityCalcData.ErrorRate = noiseRate;
        }

        if(purityCalcData.RawPurityEstimate < 0)
            purityCalcData.RawPurityEstimate = fragmentTotals.rawSampleVaf();

        purityCalcData.PurityEstimate = purityCalcData.RawPurityEstimate;

        int alleleCount = fragmentTotals.sampleAdTotal();

        double lowProbAlleleCount = calcPoissonNoiseValue(alleleCount, HIGH_PROBABILITY);

        double sampleAdjVafLow = lowProbAlleleCount / fragmentTotals.sampleDepthTotal();
        purityCalcData.PurityRangeLow = estimatedPurity(sampleAdjVafLow, noiseRate, fragmentTotals);

        double highProbAlleleCount = calcPoissonNoiseValue(alleleCount, LOW_PROBABILITY);
        double sampleAdjVafHigh = highProbAlleleCount / fragmentTotals.sampleDepthTotal();
        purityCalcData.PurityRangeHigh = estimatedPurity(sampleAdjVafHigh, noiseRate, fragmentTotals);

        if(purityCalcData.PurityRangeHigh < 0)
            purityCalcData.PurityRangeHigh = fragmentTotals.rawSampleVaf();

        if(purityCalcData.PurityRangeLow < 0)
            purityCalcData.PurityRangeLow = fragmentTotals.rawSampleVaf();

        ClonalityModel model = null;

        if(VafPeakModel.canUseModel(fragmentTotals, purityCalcData))
        {
            model = new VafPeakModel(mConfig, mResultsWriter, mSample, variants);
        }
        else
        {
            double medianVcn = medianVcn(variants);

            List<SomaticVariant> lowCountFilteredVariants = filterVariants(sampleId, fragmentTotals, variants, medianVcn);

            if(LowCountModel.canUseModel(sampleId, fragmentTotals, lowCountFilteredVariants))
            {
                model = new LowCountModel(mConfig, mResultsWriter, mSample, lowCountFilteredVariants);
            }
        }

        if(model != null)
        {
            purityCalcData.Clonality = model.calculate(sampleId, fragmentTotals, purityCalcData);

            if(isRecomputed(purityCalcData.Clonality.Method))
            {
                double lowRatio = purityCalcData.PurityEstimate > 0 ? purityCalcData.PurityRangeLow / purityCalcData.PurityEstimate : 1;
                double highRatio = purityCalcData.PurityEstimate > 0 ? purityCalcData.PurityRangeHigh / purityCalcData.PurityEstimate : 1;

                if(convertVafToPurity(purityCalcData.Clonality.Method))
                {
                    purityCalcData.PurityEstimate = estimatedPurity(purityCalcData.Clonality.Vaf, noiseRate, fragmentTotals);
                    purityCalcData.PurityRangeLow = estimatedPurity(purityCalcData.Clonality.VafLow, noiseRate, fragmentTotals);
                    purityCalcData.PurityRangeHigh = estimatedPurity(purityCalcData.Clonality.VafHigh, noiseRate, fragmentTotals);
                }
                else
                {
                    purityCalcData.PurityEstimate = purityCalcData.Clonality.Vaf;
                    purityCalcData.PurityRangeLow = purityCalcData.Clonality.VafLow;
                    purityCalcData.PurityRangeHigh = purityCalcData.Clonality.VafHigh;
                }

                purityCalcData.PurityRangeLow *= lowRatio;
                purityCalcData.PurityRangeHigh *= highRatio;
            }
        }

        // report final probability as min of Dual and Normal Prob
        double expectedDualNoiseFragments = mConfig.noiseRate(true) * sampleDualDP;
        purityCalcData.DualProbability = estimatedProbability(sampleDualAD, expectedDualNoiseFragments);

        // CT_LOGGER.info(format("patient(%s) sample(%s) sampleTotalFrags(%d) noise(%.1f) LOD(%.6f)",
        //        mSample.PatientId, sampleId, sampleDepthTotal, allFragsNoise, lodFragsResult.EstimatedPurity));

        return new SomaticPurityResult(true, totalVariantCount, chipVariants, fragmentTotals, umiTypeCounts, purityCalcData);
    }

    public double getBqrErrorRate(final SomaticVariant variant)
    {
        return mBqrAdjustment.calcErrorRate(variant.TriNucContext, variant.Alt);
    }

    private FragmentTotals calculateThresholdValues(
            final String sampleId, final PurityCalcData purityCalcData, final List<SomaticVariant> variants, int qualThreshold)
    {
        List<BqrContextData> filteredBqrData = mBqrAdjustment.getThresholdBqrData(qualThreshold);

        if(filteredBqrData.isEmpty())
            return null;

        double rawBqrErrorRate = BqrAdjustment.calcErrorRate(filteredBqrData);

        FragmentTotals fragmentTotals = new FragmentTotals();

        double depthWeightedErrorRate = 0;

        for(SomaticVariant variant : variants)
        {
            if(!hasVariantContext(filteredBqrData, variant.TriNucContext, variant.Alt))
                continue;

            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            double varBqrErrorRate = getBqrErrorRate(variant);
            sampleFragData.setBqrErrorRate(varBqrErrorRate);

            depthWeightedErrorRate += varBqrErrorRate * sampleFragData.Depth;

            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            fragmentTotals.addVariantData(
                    variant.CopyNumber, variant.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, sampleFragData.QualTotal);
        }

        if(fragmentTotals.sampleDepthTotal() == 0)
            return null;

        double bqrErrorRate = depthWeightedErrorRate / fragmentTotals.sampleDepthTotal();

        double probability = estimatedProbability(fragmentTotals, bqrErrorRate);

        double lodPurityEstimate = calcLimitOfDetection(fragmentTotals, bqrErrorRate);

        /*
        CT_LOGGER.debug(format("patient(%s) sample(%s) errorRatePM(%.1f) variants(%d) frags(dp=%d ad=%d)  probability(%.6f) LOD(%.6f)",
                mSample.PatientId, sampleId, bqrErrorRate * 1_000_000, fragmentTotals.variantCount(),
                fragmentTotals.sampleDepthTotal(), fragmentTotals.sampleAdTotal(), probability, lodPurityEstimate));
        */

        purityCalcData.BqrExtraInfo.add(format(
                "BqrThreshold(%d variants=%d AD=%d DP=%d prob=%s lod=%s)",
                qualThreshold, fragmentTotals.variantCount(), fragmentTotals.sampleAdTotal(), fragmentTotals.sampleDepthTotal(),
                formatProbabilityValue(probability), formatPurityValue(lodPurityEstimate)));

        if(purityCalcData.LodPurityEstimate == CALC_NO_SET || lodPurityEstimate < purityCalcData.LodPurityEstimate)
        {
            purityCalcData.LodPurityEstimate = lodPurityEstimate;
            purityCalcData.Probability = probability;
            purityCalcData.BqrQualThreshold = qualThreshold;
            purityCalcData.ErrorRate = bqrErrorRate;
            purityCalcData.RawBqrErrorRate = rawBqrErrorRate;
            purityCalcData.RawPurityEstimate = estimatedPurity(fragmentTotals.rawSampleVaf(), bqrErrorRate, fragmentTotals);
            return fragmentTotals;
        }

        return null;
    }

    private static double medianVcn(final List<SomaticVariant> variants)
    {
        List<Double> variantCopyNumbers = variants.stream().map(x -> x.VariantCopyNumber).collect(Collectors.toList());
        Collections.sort(variantCopyNumbers);

        int medIndex = variantCopyNumbers.size() / 2;
        return variantCopyNumbers.get(medIndex);
    }
}
