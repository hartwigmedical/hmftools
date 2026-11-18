package com.hartwig.hmftools.cobalt.norm;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_BUCKET_MIN;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.FILTER_CHECK_PERCENTILE;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.FILTER_CHECK_PERCENTILE_RATIO_THRESHOLD;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.FILTER_CHECK_PERCENTILE_STD_DEV_THRESHOLD;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.FILTER_COHORT_MIN;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.FILTER_COHORT_MIN_PERCENT;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.FILTER_MEDIAN_MAX_RATIO;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.MAPPABILITY_THRESHOLD;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.MAX_WGS_ADJUST_RATIO;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.MEDIAN_PERCENTILE;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.MIN_ADJACENT_ENRICHMENT_RATIO;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.MIN_ENRICHMENT_RATIO;
import static com.hartwig.hmftools.cobalt.norm.NormConstants.WGS_PERCENTILE_DAMPEN;
import static com.hartwig.hmftools.common.utils.Doubles.median;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class Normaliser
{
    /* Calculation steps:
    - for each sample calculate interpolated median read count per GC bucket
        - only regions meeting certain criteria are used
    - also calculate interpolated median and mean of read count across all sample buckets
    - calculate an adjusted GC ratio for each sample's region
    - calculate relative enrichment as a median from each region's sample adjusted GC ratios
        - apply a min relative enrichment threshold
    */


    public static void calcSampleAdjustedRatios(final List<String> samples, final Map<String, List<RegionData>> chrRegionData)
    {
        for(int i = 0; i < samples.size(); ++i)
        {
            NormCalcData normCalcData = calcSampleAdjustedRatios(i, chrRegionData);
            applySampleNormalisation(i, chrRegionData, normCalcData);
        }
    }

    public static NormCalcData calcSampleAdjustedRatios(final int sampleIndex, final Map<String, List<RegionData>> chrRegionData)
    {
        // calculate interpolated median read count per GC bucket across filtered regions

        List<Double> sampleReadDepths = new ArrayList<>();
        double sampleReadCountTotal = 0;

        Map<Integer, List<Double>> gcBucketReadDepths = new HashMap<>();

        for(Map.Entry<String, List<RegionData>> entry : chrRegionData.entrySet())
        {
            String chromosome = entry.getKey();

            if(!useChromosome(chromosome))
            {
                continue;
            }

            List<RegionData> regions = entry.getValue();

            for(RegionData regionData : regions)
            {
                if(regionData.mappability() < MAPPABILITY_THRESHOLD)
                    continue;

                SampleRegionData sampleRegionData = regionData.getSampleData(sampleIndex);

                if(!useSampleRegionData(sampleRegionData))
                {
                    sampleRegionData.markInvalid();
                    continue;
                }

                double readCount = sampleRegionData.ReadDepth;
                sampleReadDepths.add(readCount);
                sampleReadCountTotal += readCount;

                // take the GC bucket from the sample's region data, not the GC profile
                List<Double> bucketDepths = gcBucketReadDepths.computeIfAbsent(sampleRegionData.PanelGcBucket, k -> Lists.newArrayList());

                bucketDepths.add(readCount);
            }
        }

        if(sampleReadDepths.isEmpty())
        {
            return NormCalcData.INVALID;
        }

        double sampleMeanReadCount = sampleReadCountTotal / sampleReadDepths.size();
        double sampleMedianReadCount = median(sampleReadDepths);

        Map<Integer, Double> gcBucketMedians = new HashMap<>();

        for(Map.Entry<Integer, List<Double>> entry : gcBucketReadDepths.entrySet())
        {
            int gcBucket = entry.getKey();
            double medianReadDepth = median(entry.getValue());
            gcBucketMedians.put(gcBucket, medianReadDepth);
        }

        return new NormCalcData(sampleMeanReadCount, sampleMedianReadCount, sampleReadDepths.size(), gcBucketMedians);
    }

    public static void applySampleNormalisation(
            final int sampleIndex, final Map<String, List<RegionData>> chrRegionData, final NormCalcData normCalcData)
    {
        // apply the median values into each sample's adjusted calcs
        double sampleMedianNormalisation = normCalcData.sampleMedianNormalisation();

        for(List<RegionData> regions : chrRegionData.values())
        {
            for(RegionData regionData : regions)
            {
                SampleRegionData sampleRegionData = regionData.getSampleData(sampleIndex);

                Double gcBucketMedian = normCalcData.GcBucketMedians.get(sampleRegionData.PanelGcBucket);

                if(gcBucketMedian == null || gcBucketMedian == 0)
                {
                    sampleRegionData.markInvalid();
                    continue;
                }

                double adjustedRatio = sampleMedianNormalisation * sampleRegionData.ReadDepth / gcBucketMedian;
                sampleRegionData.setAdjustedGcRatio(adjustedRatio);
            }
        }
    }

    private static boolean useChromosome(final String chromosome)
    {
        return HumanChromosome.fromString(chromosome).isAutosome();
    }

    public static boolean useSampleRegionData(final SampleRegionData sampleRegionData)
    {
        if(sampleRegionData.ReadDepth < 0)
            return false;

        return sampleRegionData.PanelGcBucket >= GC_BUCKET_MIN && sampleRegionData.PanelGcBucket <= GC_BUCKET_MAX;
    }

    public static void calcRelativeEnrichment(
            final Map<String,List<RegionData>> chrRegionData, double minEnrichmentRatio, final WgsCopyNumberPercentiles wgsPercentiles)
    {
        for(Map.Entry<String,List<RegionData>> entry : chrRegionData.entrySet())
        {
            String chromosome = entry.getKey();
            List<RegionData> regions = entry.getValue();

            for(RegionData regionData : regions)
            {
                List<Double> sampleRelativeEnrichments = new ArrayList<>(regionData.sampleCount());

                for(SampleRegionData sampleRegionData : regionData.getSamples())
                {
                    if(!sampleRegionData.isValid())
                        continue;

                    if(sampleRegionData.WgsGcRatio < 0)
                    {
                        // -1.0 is used to mark Y regions in female samples and may also arise from masked windows in WGS data
                        continue;
                    }

                    double relativeEnrichment = sampleRegionData.adjustedGcRatio() / sampleRegionData.WgsGcRatio;
                    sampleRelativeEnrichments.add(relativeEnrichment);
                }

                Collections.sort(sampleRelativeEnrichments);

                if(filterRegion(chromosome, regionData, sampleRelativeEnrichments))
                {
                    regionData.markInvalid();
                    continue;
                }

                double percentileEnrichment = findPercentileEnrichment(sampleRelativeEnrichments, chromosome, regionData, wgsPercentiles);

                if(percentileEnrichment >= minEnrichmentRatio)
                {
                    regionData.setRelativeEnrichment(percentileEnrichment);
                }
            }

            // invalid regions with significant drop-off relative to the neighbouring ones
            //
            for(int i = 1; i < regions.size() - 1; ++i)
            {
                RegionData priorRegion = regions.get(i - 1);
                RegionData regionData = regions.get(i);
                RegionData nextRegion = regions.get(i + 1);

                if(regionData.relativeEnrichment() == 0)
                    continue;

                if(priorRegion.relativeEnrichment() > 0)
                {
                    double adjacentRatio = regionData.relativeEnrichment() / priorRegion.relativeEnrichment();

                    if(adjacentRatio < MIN_ADJACENT_ENRICHMENT_RATIO)
                    {
                        regionData.markInvalid();
                        continue;
                    }
                }

                if(nextRegion.relativeEnrichment() > 0)
                {
                    double adjacentRatio = regionData.relativeEnrichment() / nextRegion.relativeEnrichment();

                    if(adjacentRatio < MIN_ADJACENT_ENRICHMENT_RATIO)
                    {
                        regionData.markInvalid();
                    }
                }
            }
        }
    }

    private static double findPercentileEnrichment(
            final List<Double> sampleEnrichments, final String chromosome, final RegionData regionData,
            final WgsCopyNumberPercentiles wgsPercentiles)
    {
        // use a percentile derived from the WGS cohort instead of the median, which can be especially important for regions and
        // genes which are typically amplified or deleted
        double specificPercentile = wgsPercentiles.getRegionPercentile(chromosome, regionData.Position);

        double median = median(sampleEnrichments);

        if(specificPercentile == MEDIAN_PERCENTILE)
            return median;

        // dampen the percentile towards the median, eg 0.25*(1-DF)+0.5*(DF) = 0.3
        double dampenedPercentile = specificPercentile * (1 - WGS_PERCENTILE_DAMPEN) + MEDIAN_PERCENTILE * WGS_PERCENTILE_DAMPEN;

        int targetIndex = (int)floor(dampenedPercentile * sampleEnrichments.size());
        double percentileEnrichment = sampleEnrichments.get(targetIndex);

        if(median == percentileEnrichment)
            return median;

        double cappedEnrichment = min(percentileEnrichment, median * (1 + MAX_WGS_ADJUST_RATIO));
        cappedEnrichment = max(cappedEnrichment, median * (1 - MAX_WGS_ADJUST_RATIO));

        return cappedEnrichment;
    }

    @VisibleForTesting
    public static boolean filterRegion(final String chromosome, final RegionData regionData, final List<Double> sampleEnrichments)
    {
        if(sampleEnrichments.isEmpty())
            return true;

        if(HumanChromosome.fromString(chromosome).isAllosome())
            return false;

        int targetIndex = (int)floor(FILTER_CHECK_PERCENTILE * sampleEnrichments.size());
        double percentile20 = sampleEnrichments.get(targetIndex);

        // 20th percentile ratio < 0.2
        if(percentile20 >= FILTER_CHECK_PERCENTILE_RATIO_THRESHOLD)
            return false;

        double median = median(sampleEnrichments);

        List<Double> filteredEnrichments = Lists.newArrayListWithCapacity(sampleEnrichments.size());

        double filteredTotal = 0;
        for(double enrichment : sampleEnrichments)
        {
            if(enrichment > 0 && enrichment <= FILTER_MEDIAN_MAX_RATIO * median)
            {
                filteredTotal += enrichment;
                filteredEnrichments.add(enrichment);
            }
        }

        int filteredCount = filteredEnrichments.size();

        if(filteredCount == 0)
            return false;

        // Number of samples contributing to the filtered median ratio / standard deviation is at least min(10, 50% of cohort samples)
        int totalSamples = sampleEnrichments.size();

        int requiredSamples = min(FILTER_COHORT_MIN, (int)round(FILTER_COHORT_MIN_PERCENT * totalSamples));

        if(filteredCount < requiredSamples)
            return false;

        double filteredMedian = median(filteredEnrichments);
        double filteredMean = filteredTotal / filteredCount;

        // Region has filtered median ratio over our MIN_ENRICHMENT_RATIO of 0.1
        if(filteredMedian < MIN_ENRICHMENT_RATIO)
            return false;

        double varianceTotal = 0;
        for(double enrichment : filteredEnrichments)
        {
            varianceTotal += pow(enrichment - filteredMean, 2);
        }

        double stdDev = sqrt(varianceTotal / filteredCount);
        double scaledStdDev = stdDev / min(filteredMedian, 1);

        // Scaled standard deviation > 0.5
        if(scaledStdDev <= FILTER_CHECK_PERCENTILE_STD_DEV_THRESHOLD)
            return false;

        CB_LOGGER.debug(format("region(%s:%d) noise filtered: count(%d/%d) median(%.3f filtMed=%.3f 20th=%.3f) scaledSD(%.3f)",
                chromosome, regionData.Position, filteredCount, totalSamples, median, filteredMedian, percentile20, scaledStdDev));

        return true;
    }
}
