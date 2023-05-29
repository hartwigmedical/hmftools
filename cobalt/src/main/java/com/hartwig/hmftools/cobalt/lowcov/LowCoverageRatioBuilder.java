package com.hartwig.hmftools.cobalt.lowcov;

import static java.lang.String.format;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Comparators;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.cobalt.CobaltUtils;
import com.hartwig.hmftools.cobalt.ratio.RatioBuilder;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import tech.tablesaw.api.*;

public class LowCoverageRatioBuilder implements RatioBuilder
{
    private final Table mLowCoverageRatios;

    public LowCoverageRatioBuilder(final Table inputRatios, int consolidationCount,
            final ChromosomePositionCodec chromosomePosCodec)
    {
        this(inputRatios,
                Objects.requireNonNull(consolidateIntoBuckets(inputRatios, consolidationCount)),
                chromosomePosCodec);
    }

    public LowCoverageRatioBuilder(final Table inputRatios,
            final @NotNull Multimap<String, LowCovBucket> consolidateBoundaries,
            final ChromosomePositionCodec chromosomePosCodec)
    {
        CB_LOGGER.info("using {} sparse consolidated buckets, from {} input ratios",
                consolidateBoundaries.size(), inputRatios.rowCount());

        mLowCoverageRatios = CobaltUtils.createRatioTable();

        // calculate median for sample
        populateLowCoverageRatio(inputRatios, consolidateBoundaries);

        chromosomePosCodec.addEncodedChrPosColumn(mLowCoverageRatios, true);
    }

    // we use on target ratios only for now
    @Override
    public Table ratios() { return mLowCoverageRatios; }

    // we create a pan window ratio by taking the median count of super windows that combine multiple windows
    private void populateLowCoverageRatio(
        final Table rawRatios, Multimap<String, LowCovBucket> consolidateBoundaries)
    {
        for (String chromosome : rawRatios.stringColumn(CobaltColumns.CHROMOSOME).unique())
        {
            List<Double> bucketGcRatios = new ArrayList<>();

            Iterator<LowCovBucket> bucketItr = consolidateBoundaries.get(chromosome).iterator();

            if (!bucketItr.hasNext())
            {
                CB_LOGGER.fatal("low cov bucket for chromosome {} not found", chromosome);
                continue;
            }

            LowCovBucket bucket = bucketItr.next();

            Table chrRatios = rawRatios.where(rawRatios.stringColumn(CobaltColumns.CHROMOSOME).isEqualTo(chromosome));

            // make sure position is sorted
            Validate.isTrue(Comparators.isInStrictOrder(chrRatios.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).asList(),
                    Comparator.naturalOrder()));

            for (Row row : chrRatios)
            {
                if (row.getDouble(CobaltColumns.RATIO) >= 0)
                {
                    if (row.getInt(CobaltColumns.POSITION) > bucket.endPosition)
                    {
                        ReadRatio unnormalizedRatio = calcConsolidatedRatio(chromosome, bucket, bucketGcRatios);

                        if (unnormalizedRatio != null)
                        {
                            Row newRow = mLowCoverageRatios.appendRow();
                            newRow.setString(CobaltColumns.CHROMOSOME, chromosome);
                            newRow.setInt(CobaltColumns.POSITION, unnormalizedRatio.position());
                            newRow.setDouble(CobaltColumns.RATIO, unnormalizedRatio.ratio());
                        }
                        bucketGcRatios.clear();

                        if (bucketItr.hasNext())
                        {
                            // move to next bucket
                            bucket = bucketItr.next();
                        }
                        else
                        {
                            // if no more bucket we just break out
                            break;
                        }
                    }
                    bucketGcRatios.add(row.getDouble(CobaltColumns.RATIO));
                }
            }

            // add the last bucket
            ReadRatio unnormalizedRatio = calcConsolidatedRatio(chromosome, bucket, bucketGcRatios);

            if (unnormalizedRatio != null)
            {
                Row newRow = mLowCoverageRatios.appendRow();
                newRow.setString(CobaltColumns.CHROMOSOME, chromosome);
                newRow.setInt(CobaltColumns.POSITION, unnormalizedRatio.position());
                newRow.setDouble(CobaltColumns.RATIO, unnormalizedRatio.ratio());
            }
        }
    }

    @Nullable
    private static ReadRatio calcConsolidatedRatio(
            final String chromosome, LowCovBucket bucket, final List<Double> windowGcRatios)
    {
        //
        double mean = windowGcRatios.stream().mapToDouble(Double::doubleValue).sum() / windowGcRatios.size();

        CB_LOGGER.debug("consolidated window: {}:{} ({} - {}), num sub windows: {}, mean ratio: {}",
                chromosome, bucket.bucketPosition, bucket.startPosition, bucket.endPosition, windowGcRatios.size(), format("%.4f", mean));

        if (mean <= 1e-10 || Double.isNaN(mean))
        {
            return null;
        }
        return ImmutableReadRatio.builder().chromosome(chromosome).position(bucket.bucketPosition).ratio(mean).build();
    }

    @Nullable
    public static Multimap<String, LowCovBucket> calcConsolidateBuckets(final Table rawRatios,
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
    static ArrayListMultimap<String, LowCovBucket> consolidateIntoBuckets(final Table rawRatios,
            final int consolidationCount)
    {
        if (consolidationCount == 1)
            return null;

        ArrayListMultimap<String, LowCovBucket> boundaries = ArrayListMultimap.create();

        for (String chromosome : rawRatios.stringColumn(CobaltColumns.CHROMOSOME).unique())
        {
            List<Integer> nonMaskedPositions = rawRatios.where(
                    rawRatios.stringColumn(CobaltColumns.CHROMOSOME).isEqualTo(chromosome)
                    .and(rawRatios.doubleColumn(CobaltColumns.RATIO).isNonNegative()))
                    .intColumn(CobaltColumns.POSITION).asList();

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
        // make sure position is sorted
        Validate.isTrue(Comparators.isInStrictOrder(windowPositions, Comparator.naturalOrder()));

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