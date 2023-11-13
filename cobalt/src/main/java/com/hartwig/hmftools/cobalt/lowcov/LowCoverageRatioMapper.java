package com.hartwig.hmftools.cobalt.lowcov;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
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
import com.hartwig.hmftools.cobalt.ratio.RatioMapper;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import tech.tablesaw.aggregate.AggregateFunctions;
import tech.tablesaw.api.*;

public class LowCoverageRatioMapper implements RatioMapper
{
    private static final String BUCKET_ID_COLUMN = "lovCovBucketId";

    private int mConsolidationCount = 0;
    private @Nullable Multimap<String, LowCovBucket> mConsolidateBoundaries;

    private final ChromosomePositionCodec mChromosomePositionCodec;


    public LowCoverageRatioMapper(int consolidationCount, final ChromosomePositionCodec chromosomePosCodec)
    {
        mConsolidationCount = consolidationCount;
        mChromosomePositionCodec = chromosomePosCodec;
    }

    public LowCoverageRatioMapper(
            final @NotNull Multimap<String, LowCovBucket> consolidateBoundaries,
            final @NotNull ChromosomePositionCodec chromosomePosCodec)
    {
        mConsolidateBoundaries = consolidateBoundaries;
        mChromosomePositionCodec = chromosomePosCodec;
    }

    // we use on target ratios only for now
    @Override
    public Table mapRatios(final Table inputRatios)
    {
        if (mConsolidateBoundaries == null)
        {
            Validate.isTrue(mConsolidationCount > 1);
            mConsolidateBoundaries = consolidateIntoBuckets(inputRatios, mConsolidationCount);
        }

        Objects.requireNonNull(mConsolidateBoundaries);
        CB_LOGGER.info("using {} sparse consolidated buckets, from {} input ratios",
                mConsolidateBoundaries.size(), inputRatios.rowCount());

        return populateLowCoverageRatio(inputRatios, mConsolidateBoundaries);
    }

    // we create a pan window ratio by taking the mean count of super windows that combine multiple windows
    @SuppressWarnings("UnstableApiUsage")
    private Table populateLowCoverageRatio(final Table rawRatios, Multimap<String, LowCovBucket> consolidateBoundaries)
    {
        // make sure the ratios chromosome code are sorted
        Validate.isTrue(Comparators.isInStrictOrder(rawRatios.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).asList(),
                Comparator.naturalOrder()));

        String chromosome = "";
        Iterator<LowCovBucket> bucketItr = null;
        LowCovBucket bucket = null;
        int bucketId = 0;

        IntColumn bucketIdCol = IntColumn.create(BUCKET_ID_COLUMN, rawRatios.rowCount());

        // also create a table of the bucket themselves
        Table bucketTable = Table.create(
                IntColumn.create(BUCKET_ID_COLUMN),
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                BooleanColumn.create("isAutosome"));

        // first step we give each row a lowCovBucketId
        for (int i = 0; i < rawRatios.rowCount(); ++i)
        {
            Row row = rawRatios.row(i);
            long encodedChrPos = row.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS);
            String rowChr = mChromosomePositionCodec.decodeChromosome(encodedChrPos);
            int pos = mChromosomePositionCodec.decodePosition(encodedChrPos);

            if (!chromosome.equals(rowChr))
            {
                // move to next chromosome
                chromosome = rowChr;
                bucketItr = consolidateBoundaries.get(chromosome).iterator();

                if (!bucketItr.hasNext())
                {
                    CB_LOGGER.error("low cov bucket for chromosome {} not found", chromosome);
                    bucket = null;
                    assert false;
                    continue;
                }

                bucket = bucketItr.next();
                ++bucketId;
                Row bucketRow = bucketTable.appendRow();
                bucketRow.setInt(BUCKET_ID_COLUMN, bucketId);
                bucketRow.setString(CobaltColumns.CHROMOSOME, chromosome);
                bucketRow.setInt(CobaltColumns.POSITION, bucket.bucketPosition);
                bucketRow.setBoolean("isAutosome", row.getBoolean("isAutosome"));
            }

            if (bucket == null)
            {
                // no bucket for whole chromosome, or we already finished last bucket
                continue;
            }

            if (row.getDouble(CobaltColumns.RATIO) >= 0)
            {
                if (pos > bucket.endPosition)
                {
                    if (bucketItr.hasNext())
                    {
                        // move to next bucket
                        bucket = bucketItr.next();
                        ++bucketId;
                        Row bucketRow = bucketTable.appendRow();
                        bucketRow.setInt(BUCKET_ID_COLUMN, bucketId);
                        bucketRow.setString(CobaltColumns.CHROMOSOME, chromosome);
                        bucketRow.setInt(CobaltColumns.POSITION, bucket.bucketPosition);
                        bucketRow.setBoolean("isAutosome", row.getBoolean("isAutosome"));
                    }
                    else
                    {
                        // no more bucket for this chromosome. Setting bucket to null will let us skip through
                        // the rest of the chromosome
                        bucket = null;
                        continue;
                    }
                }

                bucketIdCol.set(i, bucketId);
            }
        }

        // now we assign a bucket id per row, we can do with it
        rawRatios.addColumns(bucketIdCol);

        //
        Table lovCovRatio = rawRatios.summarize(
                CobaltColumns.RATIO,
                CobaltColumns.GC_CONTENT,
                AggregateFunctions.mean).by(BUCKET_ID_COLUMN);

        CB_LOGGER.debug("low cov table: {}", lovCovRatio);

        //
        // lowCovBucketId  |   Mean [gcContent]    |  Mean [bucketEncodedChrPos]  |      Mean [ratio]
        // fix up the column names
        lovCovRatio.column("Mean [gcContent]").setName(CobaltColumns.GC_CONTENT);
        lovCovRatio.column("Mean [ratio]").setName(CobaltColumns.RATIO);

        // add is mappable column, we are not sure if this is needed yet. But if we want to pass
        // consolidated ratios to gc normalisation this is needed.
        lovCovRatio.addColumns(BooleanColumn.create(CobaltColumns.IS_MAPPABLE,
                Collections.nCopies(lovCovRatio.rowCount(), true)));

        // merge in the bucket, this is required to get the bucket position
        lovCovRatio = lovCovRatio.joinOn(BUCKET_ID_COLUMN).leftOuter(bucketTable);

        // add the encoded chromosome pos columns
        mChromosomePositionCodec.addEncodedChrPosColumn(lovCovRatio, false);

        CB_LOGGER.debug("low cov table: {}", lovCovRatio);

        return lovCovRatio;
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

            CB_LOGGER.info("chromosome: {}, lov buckets count: {}", chromosome, consolidatedBuckets.size());
        }

        return boundaries;
    }

    static int calcConsolidationCount(final double medianReadDepth)
    {
        // consolidation starts when mean read depth <= 8
        double c = 80.0 / medianReadDepth;

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
    @SuppressWarnings("UnstableApiUsage")
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