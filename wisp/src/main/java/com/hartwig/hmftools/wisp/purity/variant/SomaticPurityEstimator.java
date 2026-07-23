package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.phredQualToProbability;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.probabilityToPhredQual;
import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.HIGH_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.INDEL_ERROR_RATE;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SNV_QUAL_THRESHOLDS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SYNTHETIC_TUMOR_VAF;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.FRAGMENT_DAMPENING_FACTOR;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.MAX_PURITY_TO_CLIP;
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
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityCalcs.cappedPurity;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult.INVALID_RESULT;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class SomaticPurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;
    private final SampleData mSample;
    private final BqrAdjustment mBqrAdjustment;

    public SomaticPurityEstimator(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sampleData, final BqrAdjustment bqrAdjustment)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sampleData;
        mBqrAdjustment = bqrAdjustment;
    }

    public SomaticPurityResult calculatePurity(
            final String sampleId, final List<SomaticVariant> variants, final int totalVariantCount,
            final List<SomaticVariant> outlierVariants, boolean capMaxTumorCopyNumber)
    {
        FragmentTotals fragmentTotals = new FragmentTotals();
        FragmentTotals dualFragmentTotals = new FragmentTotals();

        int sampleDualDP = 0;
        int sampleDualAD = 0;
        boolean uncertain = false;

        UmiTypeCounts umiTypeCounts = new UmiTypeCounts();

        for(SomaticVariant variant : variants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            fragmentTotals.addVariantData(
                    variant.CopyNumber, variant.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, capMaxTumorCopyNumber);

            boolean includeDual = !sampleFragData.isDualFiltered();

            if(includeDual)
            {
                dualFragmentTotals.addVariantData(
                        variant.CopyNumber, variant.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.UmiCounts.AlleleDual,
                        tumorFragData.Depth, sampleFragData.UmiCounts.TotalDual, capMaxTumorCopyNumber);

                sampleDualDP += sampleFragData.UmiCounts.TotalDual;
                sampleDualAD += sampleFragData.UmiCounts.AlleleDual;
            }

            umiTypeCounts.add(sampleFragData.UmiCounts, includeDual);
        }

        if(fragmentTotals.sampleDepthTotal() == 0)
            return INVALID_RESULT;

        if(mConfig.hasSyntheticTumor())
        {
            fragmentTotals.setTumorVafOverride(SYNTHETIC_TUMOR_VAF);
        }

        PurityCalcData purityCalcData = new PurityCalcData();

        // firstly estimate raw purity without consideration of clonal peaks
        double noiseRate;

        // calculate a limit-of-detection (LOD), being the number of fragments that would return a 99% confidence of a tumor presence
        if(!mConfig.SkipBqr)
        {
            if(!mBqrAdjustment.hasValidData())
                return INVALID_RESULT;

            for(int threshold : SNV_QUAL_THRESHOLDS)
            {
                FragmentTotals bqrThresholdFragTotals = calculateThresholdValues(sampleId, purityCalcData, variants, threshold, capMaxTumorCopyNumber);

                if(bqrThresholdFragTotals != null)
                    fragmentTotals = bqrThresholdFragTotals;
            }

            noiseRate = purityCalcData.ErrorRate;
        }
        else
        {
            noiseRate = mConfig.noiseRate(false);
            purityCalcData.RawPurityEstimate = estimatedPurity(fragmentTotals.rawSampleVaf(), noiseRate, fragmentTotals);
            purityCalcData.Probability = estimatedProbability(fragmentTotals, noiseRate);
            purityCalcData.LodPurityEstimate = calcLimitOfDetection(fragmentTotals, noiseRate);
            purityCalcData.ErrorRate = noiseRate;
        }

        if(purityCalcData.RawPurityEstimate < 0 || purityCalcData.RawPurityEstimate > MAX_PURITY_TO_CLIP)
        {
            double recalculatedPurity = estimatedPurity(
                    fragmentTotals.rawSampleVaf() * FRAGMENT_DAMPENING_FACTOR, noiseRate, fragmentTotals);

            if(recalculatedPurity < 0 || recalculatedPurity > MAX_PURITY_TO_CLIP)
            {
                CT_LOGGER.error("Sample ({}) returned invalid purity ({})", sampleId, recalculatedPurity);
                System.exit(1);
            }
            purityCalcData.RawPurityEstimate = recalculatedPurity;
            purityCalcData.Probability = 0;
            uncertain = true;
        }
        else
        {
            purityCalcData.RawPurityEstimate = min(purityCalcData.RawPurityEstimate, 1.0);
        }

        purityCalcData.PurityEstimate = purityCalcData.RawPurityEstimate;

        int alleleCount = fragmentTotals.sampleAdTotal();

        double lowProbAlleleCount = calcPoissonNoiseValue(alleleCount, HIGH_PROBABILITY);
        double sampleAdjVafLow = lowProbAlleleCount / fragmentTotals.sampleDepthTotal();

        purityCalcData.PurityRangeLow = cappedPurity(estimatedPurity(
                sampleAdjVafLow, noiseRate, fragmentTotals), purityCalcData.RawPurityEstimate);

        double highProbAlleleCount = calcPoissonNoiseValue(alleleCount, LOW_PROBABILITY);
        double sampleAdjVafHigh = highProbAlleleCount / fragmentTotals.sampleDepthTotal();

        purityCalcData.PurityRangeHigh = cappedPurity(estimatedPurity(
                sampleAdjVafHigh, noiseRate, fragmentTotals), purityCalcData.RawPurityEstimate);

        ClonalityModel model = null;

        if(!uncertain && VafPeakModel.canUseModel(fragmentTotals, purityCalcData))
        {
            model = new VafPeakModel(mConfig, mResultsWriter, mSample, variants);
        }
        else
        {
            double medianVcn = medianVcn(variants);

            List<SomaticVariant> lowCountFilteredVariants = filterVariants(sampleId, fragmentTotals, variants, medianVcn);

            if(!uncertain && LowCountModel.canUseModel(sampleId, fragmentTotals, lowCountFilteredVariants))
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
                double INVALID_PURITY = -1;

                if(convertVafToPurity(purityCalcData.Clonality.Method))
                {
                    double purityEstimate = cappedPurity(estimatedPurity(
                            purityCalcData.Clonality.Vaf, noiseRate, fragmentTotals), INVALID_PURITY);

                    if(purityEstimate != INVALID_PURITY)
                    {
                        purityCalcData.PurityEstimate = purityEstimate;
                        purityCalcData.PurityRangeLow = cappedPurity(estimatedPurity(
                                purityCalcData.Clonality.VafLow, noiseRate, fragmentTotals), purityEstimate);

                        purityCalcData.PurityRangeHigh = cappedPurity(estimatedPurity(
                                purityCalcData.Clonality.VafHigh, noiseRate, fragmentTotals), purityEstimate);
                    }
                }
                else
                {
                    double vafPeak = cappedPurity(purityCalcData.Clonality.Vaf, INVALID_PURITY);
                    if(vafPeak != INVALID_PURITY)
                    {
                        purityCalcData.PurityEstimate = vafPeak;
                        purityCalcData.PurityRangeLow = cappedPurity(purityCalcData.Clonality.VafLow, vafPeak);
                        purityCalcData.PurityRangeHigh = cappedPurity(purityCalcData.Clonality.VafHigh, vafPeak);
                    }
                }

                purityCalcData.PurityRangeLow = min(purityCalcData.PurityRangeLow * lowRatio, 1);
                purityCalcData.PurityRangeHigh = min(purityCalcData.PurityRangeHigh * highRatio, 1);
            }
        }

        // report final probability as min of Dual and Normal Prob
        double dualNoiseRate = mConfig.noiseRate(true);
        double expectedDualNoiseFragments = dualNoiseRate * sampleDualDP;
        purityCalcData.DualProbability = estimatedProbability(sampleDualAD, expectedDualNoiseFragments);
        purityCalcData.DualLodPurityEstimate = calcLimitOfDetection(dualFragmentTotals, dualNoiseRate);

        // CT_LOGGER.info(format("patient(%s) sample(%s) sampleTotalFrags(%d) noise(%.1f) LOD(%.6f)",
        //        mSample.PatientId, sampleId, sampleDepthTotal, allFragsNoise, lodFragsResult.EstimatedPurity));

        StringJoiner sjOutlier = new StringJoiner(ITEM_DELIM);

        for(SomaticVariant outlier : outlierVariants)
        {
            GenotypeFragments sampleFragData = outlier.findGenotypeData(sampleId);
            GenotypeFragments tumorFragData = outlier.findGenotypeData(mSample.TumorId);

            FragmentTotals variantFragTotals = new FragmentTotals();

            variantFragTotals.addVariantData(
                    outlier.CopyNumber, outlier.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, false);

            double impliedTF = estimatedPurity(variantFragTotals.rawSampleVaf(), noiseRate, variantFragTotals);

            sjOutlier.add(format("%s %d/%d %.2f %s",
                    outlier, sampleFragData.AlleleCount, sampleFragData.Depth, sampleFragData.vaf(), formatPurityValue(impliedTF)));
        }

        SnvFitStatus status = uncertain ? SnvFitStatus.UNCERTAIN : SnvFitStatus.OK;

        return new SomaticPurityResult(
                status, capMaxTumorCopyNumber, totalVariantCount, sjOutlier.toString(), fragmentTotals, umiTypeCounts, purityCalcData);
    }

    private double getBqrErrorRate(final SomaticVariant variant, final GenotypeFragments sampleFragData)
    {
        return getBqrErrorRate(variant, sampleFragData, ConsensusType.NONE);
    }

    public double getBqrErrorRate(final SomaticVariant variant, final GenotypeFragments sampleFragData, final ConsensusType consensusType)
    {
        double readStrandBias = max(min(sampleFragData.readStrandBias(), 1), 0);

        double errorRateForward = mBqrAdjustment.calcErrorRate(variant.TriNucContext, variant.Alt, consensusType);
        String tncReversed = Nucleotides.reverseComplementBases(variant.TriNucContext);
        String altReversed = Nucleotides.reverseComplementBases(variant.Alt);
        double errorRateReverse = mBqrAdjustment.calcErrorRate(tncReversed, altReversed, consensusType);

        double weightedPhredQual = readStrandBias * probabilityToPhredQual(errorRateForward)
                + (1 - readStrandBias) * probabilityToPhredQual(errorRateReverse);

        return phredQualToProbability((byte)round(weightedPhredQual));
    }

    private FragmentTotals calculateThresholdValues(
            final String sampleId, final PurityCalcData purityCalcData, final List<SomaticVariant> variants, int qualThreshold,
            boolean capMaxTumorCopyNumber)
    {
        List<BqrContextData> filteredBqrData = mBqrAdjustment.getThresholdBqrData(qualThreshold);

        if(filteredBqrData.isEmpty())
            return null;

        double rawBqrErrorRate = BqrAdjustment.calcErrorRate(filteredBqrData);

        FragmentTotals fragmentTotals = new FragmentTotals();

        double depthPerCnWeightedErrorRate = 0;
        double sampleDepthPerCnTotal = 0;

        for(SomaticVariant variant : variants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            double varBqrErrorRate = 0;

            if(variant.Type == VariantType.SNP)
            {
                if(!hasVariantContext(filteredBqrData, variant.TriNucContext, variant.Alt))
                    continue;

                varBqrErrorRate = getBqrErrorRate(variant, sampleFragData);
            }
            else
            {
                varBqrErrorRate = INDEL_ERROR_RATE;
            }

            sampleFragData.setBqrErrorRate(varBqrErrorRate);

            depthPerCnWeightedErrorRate += varBqrErrorRate * sampleFragData.Depth / max(variant.CopyNumber, 1);

            sampleDepthPerCnTotal += sampleFragData.Depth / max(variant.CopyNumber, 1);

            GenotypeFragments tumorFragData = variant.findGenotypeData(mSample.TumorId);

            fragmentTotals.addVariantData(
                    variant.CopyNumber, variant.VariantCopyNumber, tumorFragData.AlleleCount, sampleFragData.AlleleCount,
                    tumorFragData.Depth, sampleFragData.Depth, capMaxTumorCopyNumber);
        }

        if(fragmentTotals.sampleDepthTotal() == 0)
            return null;

        double bqrErrorRate = depthPerCnWeightedErrorRate / sampleDepthPerCnTotal;

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
