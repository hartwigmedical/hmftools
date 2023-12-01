package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_DEPTH_PERC;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_PEAK_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.SOMATIC_PEAK_MIN_VARIANTS;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.VAF_PEAK_MODEL_MIN_AVG_DEPTH;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.SOMATICS_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.SOMATIC_PEAK_FILE_ID;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityData.NO_RESULT;
import static com.hartwig.hmftools.wisp.purity.variant.ClonalityMethod.NO_PEAK;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.kde.KernelEstimator;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.wisp.common.SampleData;
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
        // CHECK - is this still required?
        // if(estimatedResult.PurityProbability > SOMATIC_PEAK_MAX_PROBABILITY)
        //    return INVALID_RESULT;

        List<Double> variantVafRatios = Lists.newArrayList();
        // List<Double> variantVafs = Lists.newArrayList();

        double sampleWad = fragmentTotals.weightedSampleDepth();
        int depthThreshold = (int)max(VAF_PEAK_MODEL_MIN_AVG_DEPTH, sampleWad * SOMATIC_PEAK_MIN_DEPTH_PERC);

        double estimateVaf = fragmentTotals.adjSampleVaf();

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

            SampleExpVaf=SampleExpAd/SampleExpDp

            SampleExpDp=((1-SampleEstTF)*SampleWAD*2+SampleEstTF*SampleWAD*CopyNumber)/2
            SampleExpAd=round(SampleEstTF*VCN/(CopyNumber*SampleEstTF+(1-SampleEstTF)*2)*SampleExpDp,2),
                                     SampleExpVaf=round(SampleExpAd/SampleExpDp,4)
            */

            if(expectedDp <= 0)
                continue;

            double expectedVaf = expectedAd / expectedDp;

            double vafRatio = sampleFragData.vaf() / expectedVaf;

            variantVafRatios.add(vafRatio);

            // variantVafs.add(sampleFragData.vaf());

            writePeakData(mResultsWriter.getSomaticPeakWriter(), mConfig, mSample, sampleId, variant, vafRatio);
        }

        if(variantVafRatios.size() < SOMATIC_PEAK_MIN_VARIANTS)
            return NO_RESULT;

        // List<VafPeak> vafPeaks = findVafPeaks(variantVafs, estimateVaf);

        List<VafPeak> vafRatioPeaks = findVafRatioPeaks(variantVafRatios, estimateVaf);

        double sampleAdjVaf = fragmentTotals.adjSampleVaf();

        if(!vafRatioPeaks.isEmpty())
        {
            VafPeak maxPeak = vafRatioPeaks.get(vafRatioPeaks.size() - 1);
            VafPeak minPeak = vafRatioPeaks.get(0);

            CT_LOGGER.debug("sample({}) filteredVars({}) found {} somatic vaf-ratio peaks, max({})",
                    sampleId, variantVafRatios.size(), vafRatioPeaks.size(), format("%.3f", maxPeak.Peak));

            for(VafPeak vafPeak : vafRatioPeaks)
            {
                CT_LOGGER.debug("sample({}) somatic vaf-ratio peak({})", sampleId, vafPeak);
            }

            // convert VAF ratio back into VAF terms
            double maxVafPeak = maxPeak.Peak * sampleAdjVaf;
            double minVafPeak = minPeak.Peak * sampleAdjVaf;

            return new ClonalityData(ClonalityMethod.VAF_PEAK, maxVafPeak, maxVafPeak, minVafPeak, maxPeak.Count, 0);
        }

        return new ClonalityData(NO_PEAK, 0, 0, 0, 0, 0);
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

    private static final double PEAK_VAF_RATIO_BUFFER = 0.1;
    private static final double DENSITY_BANDWIDTH = 0.15; // smoothing factor
    private static final double VAF_RATIO_BUCKET = 0.1; // ratio buckets
    private static final int MAX_VAF_RATIO = 8;

    private List<VafPeak> findVafRatioPeaks(final List<Double> sampleVafRatios, double estimatedVaf)
    {
        double vafRatioTotal = 0;
        double maxVafRatio = 0;
        for(double vafRatio : sampleVafRatios)
        {
            maxVafRatio = max(maxVafRatio, vafRatio);
            vafRatioTotal += vafRatio;
        }

        if(vafRatioTotal == 0)
            return Collections.emptyList();

        double densityBandwidth = DENSITY_BANDWIDTH;

        int maxVafRatioLimit = min((int)ceil(maxVafRatio), MAX_VAF_RATIO);

        int vafRatioBuckets = (int)(maxVafRatioLimit / VAF_RATIO_BUCKET);

        double[] vafRatios = new double[vafRatioBuckets + 1];

        for(int i = 0; i <= vafRatioBuckets; ++i)
        {
            vafRatios[i] = i * VAF_RATIO_BUCKET;
        }

        // double[] vafRatios = IntStream.rangeClosed(0, vafRatioBuckets).mapToDouble(x -> x * VAF_RATIO_BUCKET / maxVafRatioLimit).toArray();
        // double[] vafs = IntStream.rangeClosed(0, maxVafLimit * vafFraction).mapToDouble(x -> x / (100d * vafFraction)).toArray();

        KernelEstimator estimator = new KernelEstimator(VAF_RATIO_BUCKET, densityBandwidth);

        sampleVafRatios.forEach(x -> estimator.addValue(x, 1.0));

        double[] densities = DoubleStream.of(vafRatios).map(estimator::getProbability).toArray();

        final List<VafPeak> peakVafs = Lists.newArrayList();

        for(int i = 1; i < densities.length - 1; i++)
        {
            double density = densities[i];

            // peak must be above the estimated VAF
            if(!Doubles.greaterThan(density, densities[i - 1]) || !Doubles.greaterThan(density, densities[i + 1]))
                continue;

            double densityVaf = vafRatios[i];

            if(densityVaf < 1.0)
                continue;

            // count up observations at this density peak
            int peakCount = 0;

            for(Double variantVar : sampleVafRatios)
            {
                if(variantVar >= densityVaf - PEAK_VAF_RATIO_BUFFER && variantVar <= densityVaf + PEAK_VAF_RATIO_BUFFER)
                    ++peakCount;
            }

            if(peakCount < SOMATIC_PEAK_MIN_PEAK_VARIANTS)
                continue;

            CT_LOGGER.debug(format("somatic peak: count(%d) vafRatio(%.3f)", peakCount, densityVaf));
            peakVafs.add(new VafPeak(densityVaf, peakCount));
        }

        Collections.sort(peakVafs);

        return peakVafs;
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
            VariantContextDecorator decorator = variant.decorator();

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

    // OLD METHOD
    private static final double PEAK_VAF_BUFFER = 0.015;

    private List<VafPeak> findVafPeaks(final List<Double> sampleVafs, double estimatedVaf)
    {
        double vafTotal = 0;
        double maxVaf = 0;
        for(double variantVaf : sampleVafs)
        {
            maxVaf = max(maxVaf, variantVaf);
            vafTotal += variantVaf;
        }

        if(vafTotal == 0)
            return Collections.emptyList();

        double avgVaf = vafTotal / sampleVafs.size();

        double densityBandwidth = max(avgVaf/8, min(avgVaf/2, 0.01));

        int maxVafLimit = min((int)round(maxVaf * 100), 99);

        // VAFs will be allocated to buckets typically of 0.002 increments, so up to 500 altogether, but capped by the max observed VAF
        int vafFraction = 5;
        double[] vafs = IntStream.rangeClosed(0, maxVafLimit * vafFraction).mapToDouble(x -> x / (100d * vafFraction)).toArray();

        KernelEstimator estimator = new KernelEstimator(0.001, densityBandwidth);

        sampleVafs.forEach(x -> estimator.addValue(x, 1.0));

        double[] densities = DoubleStream.of(vafs).map(estimator::getProbability).toArray();

        final List<VafPeak> peakVafs = Lists.newArrayList();

        for(int i = 1; i < densities.length - 1; i++)
        {
            double density = densities[i];

            // peak must be above the estimated VAF
            if(!Doubles.greaterThan(density, densities[i - 1]) || !Doubles.greaterThan(density, densities[i + 1]))
                continue;

            double densityVaf = vafs[i];

            if(densityVaf < estimatedVaf)
                continue;

            // count up observations at this density peak
            int peakCount = 0;

            for(Double variantVar : sampleVafs)
            {
                if(variantVar >= densityVaf - PEAK_VAF_BUFFER && variantVar <= densityVaf + PEAK_VAF_BUFFER)
                    ++peakCount;
            }

            if(peakCount < SOMATIC_PEAK_MIN_PEAK_VARIANTS)
                continue;

            CT_LOGGER.debug(format("somatic peak: count(%d) vaf(%.3f) densityBandwidth(%.4f)", peakCount, densityVaf, densityBandwidth));
            peakVafs.add(new VafPeak(densityVaf, peakCount));
        }

        Collections.sort(peakVafs);

        return peakVafs;
    }

}
