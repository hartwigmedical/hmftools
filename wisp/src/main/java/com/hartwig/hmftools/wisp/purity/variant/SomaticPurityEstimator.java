package com.hartwig.hmftools.wisp.purity.variant;

import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_COUNT_MODEL_MIN_FRAG_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SYNTHETIC_TUMOR_VAF;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_FRAG_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.isRecomputed;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.calcLimitOfDetection;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.estimatedProbability;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.estimatedPurity;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult.INVALID_RESULT;

import java.util.List;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class SomaticPurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final SampleData mSample;

    public SomaticPurityEstimator(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;
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
                    variant.copyNumber(), tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, sampleFragData.QualTotal);

            umiTypeCounts.add(sampleFragData.UmiCounts);
            sampleDualDP += sampleFragData.UmiCounts.TotalDual;
            sampleDualAD += sampleFragData.UmiCounts.AlleleDual;
        }

        if(fragmentTotals.sampleDepthTotal() == 0)
            return INVALID_RESULT;

        double tumorPurity = purityContext.bestFit().purity();

        if(mConfig.hasSyntheticTumor())
        {
            fragmentTotals.setTumorVafOverride(SYNTHETIC_TUMOR_VAF);
        }

        PurityCalcData purityCalcData = new PurityCalcData();

        // firstly estimate raw purity without consideration of clonal peaks
        double noiseRate = mConfig.noiseRate(false);
        purityCalcData.RawPurityEstimate = estimatedPurity(
                tumorPurity, fragmentTotals.adjTumorVaf(), fragmentTotals.adjSampleVaf(), noiseRate);

        purityCalcData.PurityEstimate = purityCalcData.RawPurityEstimate;

        double weightedAvgDepth = fragmentTotals.weightedSampleDepth();

        ClonalityModel model = null;

        if(fragmentTotals.sampleTwoPlusCount() >= SOMATIC_PEAK_MIN_FRAG_VARIANTS
        && weightedAvgDepth > SOMATIC_PEAK_MIN_AVG_DEPTH)
        {
            model = new VafPeakModel(mConfig, mResultsWriter, mSample, variants);
        }
        else if(fragmentTotals.sampleOneFragmentCount() + fragmentTotals.sampleTwoPlusCount() >= LOW_COUNT_MODEL_MIN_FRAG_VARIANTS
        && weightedAvgDepth < LOW_COUNT_MODEL_MIN_AVG_DEPTH)
        {
            model = new LowCountModel(mConfig, mResultsWriter, mSample, variants);
        }

        if(model != null)
        {
            purityCalcData.Clonality = model.calculate(sampleId, fragmentTotals, purityCalcData.RawPurityEstimate);

            if(isRecomputed(purityCalcData.Clonality.Method))
            {
                purityCalcData.PurityEstimate = estimatedPurity(
                        tumorPurity, fragmentTotals.adjTumorVaf(), purityCalcData.Clonality.Vaf, noiseRate);

                purityCalcData.PurityRangeLow = estimatedPurity(
                        tumorPurity, fragmentTotals.adjTumorVaf(), purityCalcData.Clonality.VafLow, noiseRate);

                purityCalcData.PurityRangeHigh = estimatedPurity(
                        tumorPurity, fragmentTotals.adjTumorVaf(), purityCalcData.Clonality.VafHigh, noiseRate);
            }
        }

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 95% confidence of a tumor presence
        purityCalcData.LodPurityEstimate = calcLimitOfDetection(fragmentTotals, tumorPurity, noiseRate);

        // report final probability as min of Dual and Normal Prob
        double expectedDualNoiseFragments = mConfig.noiseRate(true) * sampleDualDP;
        purityCalcData.DualProbability = estimatedProbability(sampleDualAD, expectedDualNoiseFragments);

        // CT_LOGGER.info(format("patient(%s) sample(%s) sampleTotalFrags(%d) noise(%.1f) LOD(%.6f)",
        //        mSample.PatientId, sampleId, sampleDepthTotal, allFragsNoise, lodFragsResult.EstimatedPurity));

        return new SomaticPurityResult(true, totalVariantCount, chipVariants, fragmentTotals, umiTypeCounts, purityCalcData);
    }
}
