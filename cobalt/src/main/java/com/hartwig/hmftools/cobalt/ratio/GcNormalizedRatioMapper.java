package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.Validate;

import tech.tablesaw.aggregate.AggregateFunctions;
import tech.tablesaw.aggregate.NumericAggregateFunction;
import tech.tablesaw.api.*;

public class GcNormalizedRatioMapper implements RatioMapper
{
    private static final int MIN_BUCKET = 20;
    private static final int MAX_BUCKET = 60;

    private final boolean mUseInterpolatedMedian;

    private Table mGCMedianReadCount;
    private double mSampleMedianReadCount;
    private double mSampleMeanReadCount;

    // apply gc normalisation, the input ratios must have chromosome, position, ratio, gcBucket, isMappable
    public GcNormalizedRatioMapper(boolean useInterpolatedMedian)
    {
        mUseInterpolatedMedian = useInterpolatedMedian;
    }

    @Override
    public Table mapRatios(final Table inputRatios)
    {
        CB_LOGGER.info("Applying ratio gc normalization");

        // add a gc bucket column if not already have one
        if (!inputRatios.containsColumn(CobaltColumns.GC_BUCKET))
        {
            inputRatios.addColumns(inputRatios.doubleColumn(CobaltColumns.GC_CONTENT)
                    .multiply(100).round().asIntColumn().setName(CobaltColumns.GC_BUCKET));
        }

        // create a gc normalisation df

        // skipped masked regions
        Table gcMedianCalcDf = inputRatios.where(
                inputRatios.doubleColumn(CobaltColumns.RATIO).isGreaterThan(0.0) // TODO: change to >= 0.0
                        .and(inputRatios.intColumn(CobaltColumns.GC_BUCKET).isBetweenInclusive(MIN_BUCKET, MAX_BUCKET))
                        .and(inputRatios.booleanColumn(CobaltColumns.IS_MAPPABLE).asSelection())
                        .and(inputRatios.booleanColumn(CobaltColumns.IS_AUTOSOME).asSelection()));

        NumericAggregateFunction aggFunc;

        if (mUseInterpolatedMedian)
        {
            aggFunc = new NumericAggregateFunction("Interpolated Median")
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
            aggFunc = AggregateFunctions.median;
        }

        // get the sample median and mean
        mSampleMedianReadCount = aggFunc.summarize(gcMedianCalcDf.doubleColumn(CobaltColumns.RATIO));
        mSampleMeanReadCount = gcMedianCalcDf.doubleColumn(CobaltColumns.RATIO).mean();

        // groupby gcBucket and apply median, to create a table with columns
        // gcBucket, gcMedianCount, windowCount
        gcMedianCalcDf = gcMedianCalcDf.retainColumns(CobaltColumns.GC_BUCKET, CobaltColumns.RATIO)
                .summarize(CobaltColumns.RATIO, aggFunc, AggregateFunctions.count)
                .by(CobaltColumns.GC_BUCKET);

        CB_LOGGER.trace("sample median: {}, mean: {}, gc median calc: {}",
                mSampleMedianReadCount, mSampleMeanReadCount, gcMedianCalcDf);

        gcMedianCalcDf.column(String.format("%s [%s]", aggFunc.functionName(), CobaltColumns.RATIO)).setName("gcMedianCount");
        gcMedianCalcDf.column(String.format("Count [%s]", CobaltColumns.RATIO)).setName("windowCount");

        // merge in the gc median count
        Table ratiosWithMedianCount = inputRatios
                .where(inputRatios.booleanColumn(CobaltColumns.IS_MAPPABLE).asSelection())
                .joinOn(CobaltColumns.GC_BUCKET).inner(gcMedianCalcDf);

        double medianNormalisation = mSampleMedianReadCount / mSampleMeanReadCount;

        DoubleColumn gcNormalisedRatio = ratiosWithMedianCount.doubleColumn(CobaltColumns.RATIO)
                .multiply(medianNormalisation)
                .divide(ratiosWithMedianCount.doubleColumn("gcMedianCount"))
                .map(d -> Double.isFinite(d) ? d : Double.NaN); // protect against division by 0

        ratiosWithMedianCount.replaceColumn(gcNormalisedRatio.setName(CobaltColumns.RATIO));

        // resort it, the join messes up with the ordering
        ratiosWithMedianCount = ratiosWithMedianCount.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        mGCMedianReadCount = gcMedianCalcDf;
        return ratiosWithMedianCount;
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
            medianPerBucket.put(new ImmutableGCBucket(row.getInt(CobaltColumns.GC_BUCKET)), row.getDouble("gcMedianCount"));
        }
        return new GCMedianReadCount(mSampleMeanReadCount, mSampleMedianReadCount, medianPerBucket);
    }
}
