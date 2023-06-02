package com.hartwig.hmftools.cobalt.lowcov;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.cobalt.ratio.RatioBuilder;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LowCoverageRatioBuilder implements RatioBuilder
{
    private final ArrayListMultimap<Chromosome,ReadRatio> mLowCoverageRatios;

    public LowCoverageRatioBuilder(final @NotNull Multimap<Chromosome, LowCovBucket> consolidateBoundaries, final ListMultimap<Chromosome, ReadRatio> rawRatios)
    {
        CB_LOGGER.info("using {} sparse consolidated buckets, from {} raw ratios", consolidateBoundaries.size(), rawRatios.size());

        mLowCoverageRatios = ArrayListMultimap.create();

        // calculate median for sample
        populateLowCoverageRatio(rawRatios, consolidateBoundaries);
    }

    // we use on target ratios only for now
    @Override
    public ArrayListMultimap<Chromosome, ReadRatio> ratios() { return mLowCoverageRatios; }

    // we create a pan window ratio by taking the median count of super windows that combine multiple windows
    private void populateLowCoverageRatio(
        final ListMultimap<Chromosome, ReadRatio> rawRatios, Multimap<Chromosome, LowCovBucket> consolidateBoundaries)
    {
        ArrayListMultimap<Chromosome, ReadRatio> consolidatedRatios = ArrayListMultimap.create();

        for (Chromosome chromosome : rawRatios.keySet())
        {
            List<Double> bucketGcRatios = new ArrayList<>();

            Iterator<LowCovBucket> bucketItr = consolidateBoundaries.get(chromosome).iterator();
            LowCovBucket bucket = bucketItr.next();

            // we need this to make sure we get consistent chromosome name (1 vs chr1)
            String chromosomeStr = rawRatios.get(chromosome).get(0).chromosome();

            for (ReadRatio readRatio : rawRatios.get(chromosome))
            {
                // todo: make sure this is sorted
                if (readRatio.ratio() >= 0)
                {
                    if (readRatio.position() > bucket.endPosition)
                    {
                        ReadRatio unnormalizedRatio = calcConsolidatedRatio(chromosomeStr, bucket, bucketGcRatios);

                        if (unnormalizedRatio != null)
                        {
                            consolidatedRatios.put(chromosome, unnormalizedRatio);
                        }
                        bucketGcRatios.clear();

                        if (bucketItr.hasNext())
                        {
                            // move to next bucket
                            bucket = bucketItr.next();
                        }
                        else
                        {
                            // no more bucket for this chromosome, move to next chromosome
                            break;
                        }
                    }
                    bucketGcRatios.add(readRatio.ratio());
                }
            }

            // add the last bucket
            ReadRatio unnormalizedRatio = calcConsolidatedRatio(chromosomeStr, bucket, bucketGcRatios);

            if (unnormalizedRatio != null)
            {
                consolidatedRatios.put(chromosome, unnormalizedRatio);
            }
        }

        mLowCoverageRatios.clear();
        mLowCoverageRatios.putAll(consolidatedRatios);
    }

    /*
    private void gaussianKernelRatio(
            final ListMultimap<Chromosome, ReadRatio> rawRatios, final double sigma)
    {
        ArrayListMultimap<Chromosome,ReadRatio> unnormalizedRatios = ArrayListMultimap.create();

        for (Chromosome chromosome : rawRatios.keySet())
        {
            int firstWindowStart = -1;
            int lastWindowStart = -1;

            List<Double> bucketGcRatios = new ArrayList<>();

            // we need this to make sure we get consistent chromosome name (1 vs chr1)
            String chromosomeStr = rawRatios.get(chromosome).get(0).chromosome();

            // TODO remove outlier

            for (ReadRatio readRatio : rawRatios.get(chromosome))
            {
                // todo: make sure this is sorted
                if (readRatio.ratio() >= 0)
                {
                    if (firstWindowStart == -1)
                    {
                        firstWindowStart = readRatio.position();
                    }

                    lastWindowStart = readRatio.position();
                    bucketGcRatios.add(readRatio.ratio());

                    if (bucketGcRatios.size() == consolidationSize)
                    {
                        if (firstWindowStart != -1)
                        {
                            ReadRatio unnormalizedRatio = calcConsolidatedRatio(
                                    chromosomeStr, firstWindowStart, lastWindowStart + CobaltConstants.WINDOW_SIZE, bucketGcRatios);

                            if (unnormalizedRatio != null)
                            {
                                unnormalizedRatios.put(chromosome, unnormalizedRatio);
                            }
                        }

                        firstWindowStart = -1;
                        bucketGcRatios.clear();
                    }
                }
            }

            ReadRatio unnormalizedRatio = calcConsolidatedRatio(
                    chromosomeStr, firstWindowStart, lastWindowStart + CobaltConstants.WINDOW_SIZE, bucketGcRatios);

            if (unnormalizedRatio != null)
                unnormalizedRatios.put(chromosome, unnormalizedRatio);
        }

        // now we want to normalise all of those off target gc ratios by the median
        List<Double> values = new ArrayList<>();
        unnormalizedRatios.values().forEach(x -> values.add(x.ratio()));
        double median = Doubles.median(values);

        CB_LOGGER.debug("normalizing {} off target windows ratio by median: {}", unnormalizedRatios.size(), median);

        mLowCoverageRatios.clear();

        //for ((key, value) in unnormalizedRatios.entries())
        for (Map.Entry<Chromosome,ReadRatio> entry : unnormalizedRatios.entries())
        {
            ReadRatio readRatio = entry.getValue();
            double normalizedRatio = readRatio.ratio() / median;

            mLowCoverageRatios.put(entry.getKey(), ImmutableReadRatio.builder().from(readRatio).ratio(normalizedRatio).build());
        }
    }*/

    @Nullable
    private static ReadRatio calcConsolidatedRatio(
            final String chromosome, LowCovBucket bucket, final List<Double> windowGcRatios)
    {
        //
        double mean = windowGcRatios.stream().mapToDouble(Double::doubleValue).sum() / windowGcRatios.size();

        CB_LOGGER.debug("consolidated window: {}:{} ({} - {}), num sub windows: {}, mean ratio: {}",
                chromosome, bucket.bucketPosition, bucket.startPosition, bucket.endPosition, windowGcRatios.size(), format("%.4f", mean));

        return !Double.isNaN(mean) ? ImmutableReadRatio.builder().chromosome(chromosome).position(bucket.bucketPosition).ratio(mean).build() : null;
    }

    @Nullable
    public static Multimap<Chromosome, LowCovBucket> calcConsolidateBuckets(final ListMultimap<Chromosome, ReadRatio> rawRatios,
            final double medianReadCount)
    {
        int consolidationCount = calcConsolidationCount(medianReadCount);

        if (consolidationCount == 1)
        {
            CB_LOGGER.info("median read count: {}, not using sparse consolidation", medianReadCount);
            return null;
        }

        CB_LOGGER.info("median read count: {}, sparse consolidation count: {}",
                medianReadCount, consolidationCount);

        return consolidateIntoBuckets(rawRatios, consolidationCount);
    }

    // given the consolidation count, which is the number of 1k window we want in each bucket, we go through the windows and
    // and find the ranges of the consolidated buckets. We do this to skip through windows with invalid ratios.
    @Nullable
    static ArrayListMultimap<Chromosome, LowCovBucket> consolidateIntoBuckets(final ListMultimap<Chromosome, ReadRatio> rawRatios,
            final int consolidationCount)
    {
        if (consolidationCount == 1)
            return null;

        ArrayListMultimap<Chromosome, LowCovBucket> boundaries = ArrayListMultimap.create();

        for (Chromosome chromosome : rawRatios.keySet())
        {
            List<Integer> nonMaskedPositions = rawRatios.get(chromosome).stream()
                    .filter(o -> o.ratio() >= 0.0)
                    .mapToInt(ReadRatio::position).boxed()
                    .collect(Collectors.toList());

            List<LowCovBucket> consolidatedBuckets = consolidateIntoBuckets(nonMaskedPositions, consolidationCount);

            boundaries.putAll(chromosome, consolidatedBuckets);
        }

        return boundaries;
    }

    static int calcConsolidationCount(final double medianReadCount)
    {
        // consolidation starts when mean read count <= 50 reads
        double c = 500.0 / medianReadCount;

        if (c < 10.0)
        {
            return 1;
        }

        // max 1000
        if (c >= 1_000)
        {
            return 1_000;
        }

        // round to one significant digit
        double roundBy = Math.pow(10, Math.floor(Math.log10(c)));
        return (int) (Math.round(c / roundBy) * roundBy);
    }

    private static int roundDownToWindowBoundary(double p)
    {
        return (int) (Math.floor(p / CobaltConstants.WINDOW_SIZE) * CobaltConstants.WINDOW_SIZE) + 1;
    }

    // given the list of non masked windows, get the list of consolidated buckets
    static List<LowCovBucket> consolidateIntoBuckets(List<Integer> windowPositions, int consolidationCount)
    {
        List<LowCovBucket> buckets = new ArrayList<>();

        if (windowPositions.isEmpty())
        {
            return buckets;
        }

        int windowCount = 0;
        int bucketStart = windowPositions.get(0);

        for (int i = 0; i < windowPositions.size(); ++i)
        {
            int position = windowPositions.get(i);

            if ((position - bucketStart) >= CobaltConstants.MAX_SPARSE_CONSOLIDATE_DISTANCE)
            {
                // do not let the bucket to consolidate over 3M bases. This is done to avoid consolidating
                // over centromere
                // use the last bucket
                if (i > 0)
                {
                    int lastPosition = windowPositions.get(i - 1);
                    int bucketEnd = lastPosition + CobaltConstants.WINDOW_SIZE;
                    int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
                    buckets.add(new LowCovBucket(bucketStart, bucketEnd, bucketPos));
                }

                // reset bucket start position, we do not want to put the start in the middle since
                // it would sit inside the centromere
                bucketStart = position;

                // also reset window count
                windowCount = 0;
            }

            if (windowCount == consolidationCount)
            {
                // we want to put the bucket boundary in the middle of the two windows
                int lastPosition = windowPositions.get(i - 1);
                int bucketEnd = roundDownToWindowBoundary((lastPosition + position) * 0.5);

                // bucket position is at the middle
                int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);

                buckets.add(new LowCovBucket(bucketStart, bucketEnd, bucketPos));

                // next bucket starts right after this
                bucketStart = bucketEnd + CobaltConstants.WINDOW_SIZE;

                // also reset window count
                windowCount = 0;
            }

            windowCount++;
        }

        // add a final window
        if (windowCount > 0)
        {
            int bucketEnd = windowPositions.get(windowPositions.size() - 1) + CobaltConstants.WINDOW_SIZE;

            // bucket position is at the middle
            int bucketPos = roundDownToWindowBoundary((bucketStart + bucketEnd) * 0.5);
            buckets.add(new LowCovBucket(bucketStart, bucketEnd, bucketPos));
        }

        return buckets;
    }
}