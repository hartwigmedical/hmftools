package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import static tech.tablesaw.aggregate.AggregateFunctions.count;
import static tech.tablesaw.aggregate.AggregateFunctions.median;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;
import com.hartwig.hmftools.common.utils.Doubles;

import tech.tablesaw.aggregate.NumericAggregateFunction;
import tech.tablesaw.api.*;

public class GcNormalizedRatioBuilder
{
    private static final int MIN_BUCKET = 20;
    private static final int MAX_BUCKET = 60;

    private final Table mGCMedianReadCount;
    private final double mSampleMedianReadCount;
    private final double mSampleMeanReadCount;
    private final Table mGcRatios;

    // apply gc normalisation, the input ratios must have chromosome, position, ratio, gcBucket, isMappable
    public GcNormalizedRatioBuilder(Table inputRatios, boolean useInterpolatedMedian)
    {
        CB_LOGGER.info("Applying ratio gc normalization");

        // create a gc normalisation df

        // skipped masked regions
        Table gcMedianCalcDf = inputRatios.where(
                inputRatios.doubleColumn("ratio").isNonNegative()
                .and(inputRatios.intColumn("gcBucket").isBetweenExclusive(MIN_BUCKET, MAX_BUCKET))
                .and(inputRatios.booleanColumn("isMappable").asSelection())
                .and(inputRatios.booleanColumn("isAutosome").asSelection()));

        NumericAggregateFunction aggFunc;

        if (useInterpolatedMedian)
        {
            aggFunc = new NumericAggregateFunction("interpolatedMedian")
            {
                @Override
                public Double summarize(NumericColumn<?> column)
                {
                    return Doubles.interpolatedMedian(column.asDoubleColumn().removeMissing().asList());
                }
            };
        }
        else
        {
            // use normal median
            aggFunc = median;
        }

        // get the sample median and mean
        mSampleMedianReadCount = aggFunc.summarize(gcMedianCalcDf.doubleColumn("ratio"));
        mSampleMeanReadCount = gcMedianCalcDf.doubleColumn("ratio").mean();

        // groupby gcBucket and apply median, to create a table with columns
        // gcBucket, gcMedianCount, windowCount
        gcMedianCalcDf = gcMedianCalcDf.selectColumns("gcBucket", "ratio")
                .summarize("ratio", aggFunc, count)
                .by("gcBucket");
        gcMedianCalcDf.column(1).setName("gcMedianCount");
        gcMedianCalcDf.column(2).setName("windowCount");

        // merge in the gc median count
        Table ratiosWithMedianCount = inputRatios.joinOn("gcBucket").leftOuter(gcMedianCalcDf);

        double medianNormalisation = mSampleMedianReadCount / mSampleMeanReadCount;

        DoubleColumn gcNormalisedRatio = ratiosWithMedianCount.doubleColumn(CobaltColumns.RATIO)
                .multiply(medianNormalisation)
                .divide(ratiosWithMedianCount.doubleColumn("gcMedianCount"));

        ratiosWithMedianCount.replaceColumn(gcNormalisedRatio.setName(CobaltColumns.RATIO));

        // resort it, the join messes up with the ordering
        ratiosWithMedianCount = ratiosWithMedianCount.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        mGCMedianReadCount = gcMedianCalcDf;
        mGcRatios = ratiosWithMedianCount;
    }

    public Table ratios()
    {
        return mGcRatios;
    }

    public Table gcMedianReadCountTable()
    {
        return mGCMedianReadCount;
    }

    public double getSampleMedianReadCount()
    {
        return mSampleMedianReadCount;
    }

    public double getSampleMeanReadCount()
    {
        return mSampleMeanReadCount;
    }

    // convert the gc median read count table to the object representation
    public GCMedianReadCount gcMedianReadCount()
    {
        final Map<GCBucket, Double> medianPerBucket = new HashMap<>();
        for (Row row : mGCMedianReadCount)
        {
            medianPerBucket.put(new ImmutableGCBucket(row.getInt("gcBucket")), row.getDouble("gcMedianCount"));
        }
        return new GCMedianReadCount(mSampleMeanReadCount, mSampleMedianReadCount, medianPerBucket);
    }
}
