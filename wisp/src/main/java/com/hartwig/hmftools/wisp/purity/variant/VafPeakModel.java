package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.stats.PoissonCalcs.calcPoissonNoiseValue;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.HIGH_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_BANDWIDTH_MAX;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_BANDWIDTH_MIN;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_DEPTH_PERC;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_PEAK_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_PEAK_VARIANTS_PERC;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_NTH_RATIO;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_NTH_RATIO_MIN;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.SOMATIC_PEAK_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityData.NO_RESULT;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.NO_PEAK;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.VAF_PEAK;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.DoubleStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.wisp.purity.SampleData;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class VafPeakModel extends ClonalityModel
{
    public VafPeakModel(
            final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample, final List<SomaticVariant> variants)
    {
        super(config, resultsWriter, sample,  variants);
    }

    @Override
    public ClonalityData calculate(final String sampleId, final FragmentTotals fragmentTotals, final double rawEstimatedPurity)
    {
        List<Double> variantVafRatios = Lists.newArrayList();

        double sampleWad = fragmentTotals.weightedSampleDepth();
        int depthThreshold = (int)max(SOMATIC_PEAK_MIN_AVG_DEPTH, sampleWad * SOMATIC_PEAK_MIN_DEPTH_PERC);

        for(SomaticVariant variant : mVariants)
        {
            GenotypeFragments sampleFragData = variant.findGenotypeData(sampleId);

            if(sampleFragData == null)
                continue;

            if(!canUseVariant(variant, sampleFragData, depthThreshold))
                continue;

            double copyNumber = variant.copyNumber();
            double variantCopyNumber = variant.variantCopyNumber();

            double expectedDp = ((1 - rawEstimatedPurity) * sampleWad * 2 + rawEstimatedPurity * sampleWad * copyNumber) / 2;

            double expectedAd = rawEstimatedPurity * variantCopyNumber / (copyNumber * rawEstimatedPurity + (1 - rawEstimatedPurity) * 2) * expectedDp;

            /*
            SampleExpDp = ((1 - SampleEstTF) * SampleWAD * 2 + SampleEstTF * SampleWAD * CopyNumber) / 2
            SampleExpAd = SampleEstTF * VCN / (CopyNumber * SampleEstTF + (1 - SampleEstTF) * 2) * SampleExpDp
            SampleExpVaf= SampleExpAd / SampleExpDp
            */

            if(expectedDp <= 0)
                continue;

            double expectedVaf = expectedAd / expectedDp;

            double vafRatio = sampleFragData.vaf() / expectedVaf;

            variantVafRatios.add(vafRatio);

            writePeakData(mResultsWriter.getSomaticPeakWriter(), mConfig, mSample, sampleId, variant, vafRatio);
        }

        int variantCount = variantVafRatios.size();

        if(variantCount < SOMATIC_PEAK_MIN_VARIANTS)
            return NO_RESULT;

        Collections.sort(variantVafRatios);

        double sampleAdjVaf = fragmentTotals.adjSampleVaf();

        double densityBandwidth = calculateDensityBandwidth(sampleAdjVaf, sampleWad, variantVafRatios, 1);

        VafPeak vafRatioPeak = findMaxVafRatioPeak(variantVafRatios, densityBandwidth);

        ClonalityMethod method = NO_PEAK;
        double mainPeak = 1;

        if(vafRatioPeak != null && vafRatioPeak.Peak > 1)
        {
            mainPeak = vafRatioPeak.Peak;
            method = VAF_PEAK;
        }

        int alleleCount = fragmentTotals.sampleAdTotal();

        double lowProbAlleleCount = calcPoissonNoiseValue(alleleCount, HIGH_PROBABILITY);
        double sampleAdjVafLow = fragmentTotals.adjSampleVaf(lowProbAlleleCount - alleleCount);
        double densityBandwidthLow = calculateDensityBandwidth(sampleAdjVafLow, sampleWad, variantVafRatios, 2);
        VafPeak vafRatioPeakLow = findMaxVafRatioPeak(variantVafRatios, densityBandwidthLow);

        double highProbAlleleCount = calcPoissonNoiseValue(alleleCount, LOW_PROBABILITY);
        double sampleAdjVafHigh = fragmentTotals.adjSampleVaf(highProbAlleleCount - alleleCount);
        double densityBandwidthHigh = calculateDensityBandwidth(sampleAdjVafHigh, sampleWad, variantVafRatios, 0.5);
        VafPeak vafRatioPeakHigh = findMaxVafRatioPeak(variantVafRatios, densityBandwidthHigh);

        double peakHigh = 0;
        double peakLow = 0;

        if(vafRatioPeakLow != null && vafRatioPeakHigh != null)
        {
            peakHigh = max(max(vafRatioPeakLow.Peak, vafRatioPeakHigh.Peak), mainPeak);
            peakLow = min(min(vafRatioPeakLow.Peak, vafRatioPeakHigh.Peak), mainPeak);
        }
        else if(vafRatioPeakLow != null)
        {
            peakHigh = max(vafRatioPeakLow.Peak, mainPeak);
            peakLow = min(vafRatioPeakLow.Peak, mainPeak);
        }
        else if(vafRatioPeakHigh != null)
        {
            peakHigh = max(vafRatioPeakHigh.Peak, mainPeak);
            peakLow = min(vafRatioPeakHigh.Peak, mainPeak);
        }

        return new ClonalityData(
                method,
                mainPeak * sampleAdjVaf,
                peakLow * sampleAdjVafLow,
                peakHigh * sampleAdjVafHigh,
                vafRatioPeak != null ? vafRatioPeak.Count : 0,
                0,
                densityBandwidth, densityBandwidthLow, densityBandwidthHigh);
    }

    private double calculateDensityBandwidth(
            final double sampleAdjVaf, final double weightedSampleDepth, final List<Double> sampleVafRatios, final double multiplier)
    {
        /*
        nthRatio = round(pmax(3,0.08*nrow(peakVars)))
        calcBW = peakVars %>% arrange(-VafRatio) %>% summarise(NthRatio=nth(VafRatio-1,nthRatio))
        minBW=10/(adjSampleVaf*weightedAvgDepth)
        finalBW=pmin(pmax(0.2,pmax(calcBW$NthRatio/4,minBW)),3)
         */

        int variantCount = sampleVafRatios.size();
        int nthItem = (int)round(max(SOMATIC_PEAK_NTH_RATIO_MIN, SOMATIC_PEAK_NTH_RATIO * variantCount));

        double nthRatio = sampleVafRatios.get(variantCount - nthItem);
        double minBandwidth = 10 / (sampleAdjVaf * weightedSampleDepth);

        double densityBandwidth = max((nthRatio - 1) / 4, minBandwidth);

        densityBandwidth *= multiplier;

        return min(max(SOMATIC_PEAK_BANDWIDTH_MIN, densityBandwidth), SOMATIC_PEAK_BANDWIDTH_MAX);
    }

    private static final double PEAK_VAF_RATIO_BUFFER = 0.1;
    private static final double VAF_RATIO_BUCKET = 0.1; // ratio buckets
    private static final int MAX_VAF_RATIO = 20;

    private VafPeak findMaxVafRatioPeak(final List<Double> sampleVafRatios, double densityBandwidth)
    {
        double maxVafRatio = 0;
        for(double vafRatio : sampleVafRatios)
        {
            maxVafRatio = max(maxVafRatio, vafRatio);
        }

        int maxVafRatioLimit = min((int)ceil(maxVafRatio), MAX_VAF_RATIO);

        int vafRatioBuckets = (int)(maxVafRatioLimit / VAF_RATIO_BUCKET);

        double[] vafRatios = new double[vafRatioBuckets + 1];

        for(int i = 0; i <= vafRatioBuckets; ++i)
        {
            vafRatios[i] = i * VAF_RATIO_BUCKET;
        }

        KernelEstimator estimator = new KernelEstimator(VAF_RATIO_BUCKET, densityBandwidth);

        sampleVafRatios.forEach(x -> estimator.addValue(x, 1.0));

        double[] densities = DoubleStream.of(vafRatios).map(estimator::getProbability).toArray();

        double densityTotal = Arrays.stream(densities).sum();
        int ratioCount = sampleVafRatios.size();

        final List<VafPeak> vafPeaks = Lists.newArrayList();

        int peakStart = 1;
        int lastPeakIndex = -1;

        for(int i = 1; i < densities.length - 1; i++)
        {
            double density = densities[i];

            // identify a trough
            if(Doubles.lessThan(density, densities[i - 1]) && Doubles.lessThan(density, densities[i + 1]))
            {
                peakStart = i;
                continue;
            }

            // identify a peak
            if(!Doubles.greaterThan(density, densities[i - 1]) || !Doubles.greaterThan(density, densities[i + 1]))
                continue;

            double vafRatio = vafRatios[i];

            if(peakStart < 0)
                continue;

            if(lastPeakIndex > 0 && peakStart < lastPeakIndex)
                continue;

            // sum density observations at this density peak
            double peakDensityTotal = 0;
            int peakEnd = -1;

            for(int j = peakStart; j < densities.length; ++j)
            {
                double varDensity = densities[j];

                // break if peak ends (ie an up-tick)
                peakEnd = j;
                if(j > i && j < densities.length - 1 && Doubles.lessThan(varDensity, densities[j + 1]))
                {
                    break;
                }

                peakDensityTotal += varDensity;
            }

            double impliedPeakVarCount = peakDensityTotal / densityTotal * ratioCount;
            double impliedPeakVarPerc = impliedPeakVarCount / ratioCount;

            if(impliedPeakVarCount < SOMATIC_PEAK_MIN_PEAK_VARIANTS || impliedPeakVarPerc < SOMATIC_PEAK_MIN_PEAK_VARIANTS_PERC)
                continue;

            int peakCount = (int)round(impliedPeakVarCount);

            /*
            int peakCount = 0;
            for(Double variantVar : sampleVafRatios)
            {
                if(variantVar >= vafRatio - PEAK_VAF_RATIO_BUFFER && variantVar <= vafRatio + PEAK_VAF_RATIO_BUFFER)
                    ++peakCount;
            }

            if(peakCount < SOMATIC_PEAK_MIN_PEAK_VARIANTS)
                continue;
            */

            CT_LOGGER.trace(format("somatic vafRatio peak(%.3f @ %d) count(%d) range(%d - %d)",
                    vafRatio, i, peakCount, peakStart, peakEnd));

            vafPeaks.add(new VafPeak(vafRatio, peakCount));
            lastPeakIndex = peakStart;
            peakStart = -1;
        }

        if(vafPeaks.isEmpty())
            return null;

        // return the highest
        return vafPeaks.get(vafPeaks.size() - 1);
    }

    private class VafPeak implements Comparable<VafPeak>
    {
        public final double Peak;
        public final int Count;

        public VafPeak(final double peak, final int count)
        {
            Peak = peak;
            Count = count;
        }

        @Override
        public int compareTo(final VafPeak other)
        {
            if(Peak == other.Peak)
                return 0;

            return Peak < other.Peak ? -1 : 1;
        }

        public String toString() { return format("%.3f=%d", Peak, Count); }
    }

    public static BufferedWriter initialiseSomaticPeakWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(SOMATIC_PEAK_FILE_ID);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add("VariantInfo").add("VafRatio");
            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    private static synchronized void writePeakData(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final SomaticVariant variant, final double vafRatio)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonFields(sj, config, sampleData, sampleId);

            sj.add(format("%s %d %s>%s", variant.Chromosome, variant.Position, variant.Ref, variant.Alt));
            sj.add(format("%.4f", vafRatio));

            writer.write(sj.toString());

            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
            System.exit(1);
        }
    }

    private boolean canUseVariant(final SomaticVariant variant, final GenotypeFragments sampleFragData, int depthThreshold)
    {
        return useVariant(variant, sampleFragData)
            && sampleFragData.UmiCounts.totalCount() >= depthThreshold
            && sampleFragData.UmiCounts.alleleCount() >= 1;
    }
}
