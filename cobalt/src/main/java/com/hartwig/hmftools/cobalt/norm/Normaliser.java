package com.hartwig.hmftools.cobalt.utils;

import static com.hartwig.hmftools.cobalt.utils.NormConstants.GC_BUCKET_MAX;
import static com.hartwig.hmftools.cobalt.utils.NormConstants.GC_BUCKET_MIN;
import static com.hartwig.hmftools.cobalt.utils.NormConstants.MAPPABILITY_THRESHOLD;
import static com.hartwig.hmftools.common.utils.Doubles.interpolatedMedian;
import static com.hartwig.hmftools.common.utils.Doubles.median;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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

    public static void calcSampleAdjustedRatios(final List<String> samples, final Map<String,List<RegionData>> chrRegionData)
    {
        for(int i = 0; i < samples.size(); ++i)
        {
            calcSampleAdjustedRatios(i, chrRegionData);
        }
    }

    private static void calcSampleAdjustedRatios(final int sampleIndex, final Map<String,List<RegionData>> chrRegionData)
    {
        // calculate interpolated median read count per GC bucket across filtered regions

        // gather up the sample's region data entries

        List<Double> sampleReadCounts = Lists.newArrayList();
        double sampleReadCountTotal = 0;

        Map<Integer, List<Double>> gcBucketReadCounts = Maps.newHashMap();

        for(Map.Entry<String, List<RegionData>> entry : chrRegionData.entrySet())
        {
            String chromosome = entry.getKey();

            if(!useChromosome(chromosome))
                continue;

            List<RegionData> regions = entry.getValue();

            for(RegionData regionData : regions)
            {
                SampleRegionData sampleRegionData = regionData.getSampleData(sampleIndex);

                if(!useRegion(regionData, sampleRegionData))
                    continue;

                double readCount = sampleRegionData.ReadCount;
                sampleReadCounts.add(readCount);
                sampleReadCountTotal += readCount;

                List<Double> bucketCounts = gcBucketReadCounts.get(regionData.gcBucket());

                if(bucketCounts == null)
                {
                    bucketCounts = Lists.newArrayList();
                    gcBucketReadCounts.put(regionData.gcBucket(), bucketCounts);
                }

                bucketCounts.add(readCount);
            }
        }

        if(sampleReadCounts.isEmpty())
            return;

        double sampleMeanReadCount = sampleReadCountTotal / sampleReadCounts.size();
        double sampleMedianReadCount = interpolatedMedian(sampleReadCounts);
        double sampleMedianNormalisation = sampleMedianReadCount / sampleMeanReadCount;

        Map<Integer,Double> gcBucketMedians = Maps.newHashMap();

        for(Map.Entry<Integer, List<Double>> entry : gcBucketReadCounts.entrySet())
        {
            int gcBucket = entry.getKey();
            double medianReadCount = interpolatedMedian(entry.getValue());
            gcBucketMedians.put(gcBucket, medianReadCount);
        }

        // now combine these median values into each sample's adjusted calcs
        for(List<RegionData> regions : chrRegionData.values())
        {
            for(RegionData regionData : regions)
            {
                SampleRegionData sampleRegionData = regionData.getSampleData(sampleIndex);

                Double gcBucketMedian = gcBucketMedians.get(regionData.gcBucket());

                if(gcBucketMedian == null || gcBucketMedian == 0)
                    continue;

                double adjustedRatio = sampleMedianNormalisation * sampleRegionData.ReadCount / gcBucketMedian;

                sampleRegionData.setAdjustedGcRatio(adjustedRatio);
            }
        }
    }

    private static boolean useChromosome(final String chromosome)
    {
        return HumanChromosome.fromString(chromosome).isAutosome();
    }

    public static boolean useRegion(final RegionData regionData, final SampleRegionData sampleRegionData)
    {
        if(regionData.mappability() < MAPPABILITY_THRESHOLD)
            return false;

        if(sampleRegionData.ReadCount < 0)
            return false;

        if(regionData.gcBucket() < GC_BUCKET_MIN || regionData.gcBucket() > GC_BUCKET_MAX)
            return false;

        return true;
    }

    public static void calcRelativeEnrichment(final Map<String,List<RegionData>> chrRegionData, double minEnrichmentRatio)
    {
        for(List<RegionData> regions : chrRegionData.values())
        {
            for(RegionData regionData : regions)
            {
                List<Double> sampleRelativeEnrichment = Lists.newArrayListWithCapacity(regionData.sampleCount());

                for(SampleRegionData sampleRegionData : regionData.getSamples())
                {
                    double relativeEnrichment = sampleRegionData.GcRatioWgs > 0 ?
                            sampleRegionData.adjustedGcRatio() / sampleRegionData.GcRatioWgs : 0;

                    sampleRelativeEnrichment.add(relativeEnrichment);
                }

                double medianEnrichment = median(sampleRelativeEnrichment);

                if(medianEnrichment >= minEnrichmentRatio)
                    regionData.setRelativeEnrichment(medianEnrichment);
            }
        }
    }
}
