package com.hartwig.hmftools.cobalt.ratio;

import static java.lang.Math.round;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_RATIO_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.GC_RATIO_MIN;

import java.util.HashMap;
import java.util.Map;

import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.gc.GCBucket;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadDepth;
import com.hartwig.hmftools.common.genome.gc.ImmutableGCBucket;

import tech.tablesaw.aggregate.AggregateFunctions;
import tech.tablesaw.aggregate.NumericAggregateFunction;
import tech.tablesaw.api.*;

public class GcNormalizedRatioMapper implements RatioMapper
{
    private Table mGCMedianReadDepth;
    private double mSampleMedianReadDepth;
    private double mSampleMeanReadDepth;

    // apply gc normalisation, the input ratios must have chromosome, position, ratio, gcBucket, isMappable
    public GcNormalizedRatioMapper() {}

    @Override
    public Table mapRatios(final Table inputRatios)
    {
        CB_LOGGER.info("applying ratio GC normalisation");

        // add a gc bucket column if not already have one
        if (!inputRatios.containsColumn(CobaltColumns.GC_BUCKET))
        {
            inputRatios.addColumns(inputRatios.doubleColumn(CobaltColumns.GC_CONTENT)
                    .multiply(100).round().asIntColumn().setName(CobaltColumns.GC_BUCKET));
        }

        // create a gc normalisation df

        int gcRatioBucketMin = (int)round(GC_RATIO_MIN * 100);
        int gcRatioBucketMax = (int)round(GC_RATIO_MAX * 100);

        // skipped masked regions
        Table gcMedianCalcDf = inputRatios.where(
                inputRatios.doubleColumn(CobaltColumns.RATIO).isGreaterThan(0.0) // TODO: change to >= 0.0
                        .and(inputRatios.intColumn(CobaltColumns.GC_BUCKET).isBetweenInclusive(gcRatioBucketMin, gcRatioBucketMax))
                        .and(inputRatios.booleanColumn(CobaltColumns.IS_MAPPABLE).asSelection())
                        .and(inputRatios.booleanColumn(CobaltColumns.IS_AUTOSOME).asSelection()));

        NumericAggregateFunction aggFunc = AggregateFunctions.median;

        // get the sample median and mean
        mSampleMedianReadDepth = aggFunc.summarize(gcMedianCalcDf.doubleColumn(CobaltColumns.RATIO));
        mSampleMeanReadDepth = gcMedianCalcDf.doubleColumn(CobaltColumns.RATIO).mean();

        // groupby gcBucket and apply median, to create a table with columns
        // gcBucket, gcMedianCount, windowCount
        gcMedianCalcDf = gcMedianCalcDf.retainColumns(CobaltColumns.GC_BUCKET, CobaltColumns.RATIO)
                .summarize(CobaltColumns.RATIO, aggFunc, AggregateFunctions.count)
                .by(CobaltColumns.GC_BUCKET);

        CB_LOGGER.trace("sample median: {}, mean: {}, gc median calc: {}", mSampleMedianReadDepth, mSampleMeanReadDepth, gcMedianCalcDf);

        gcMedianCalcDf.column(String.format("Median [%s]", CobaltColumns.RATIO)).setName("gcMedianCount");
        gcMedianCalcDf.column(String.format("Count [%s]", CobaltColumns.RATIO)).setName("windowCount");

        // merge in the gc median count
        Table ratiosWithMedianCount = inputRatios
                .where(inputRatios.booleanColumn(CobaltColumns.IS_MAPPABLE).asSelection())
                .joinOn(CobaltColumns.GC_BUCKET).inner(gcMedianCalcDf);

        double medianNormalisation = mSampleMedianReadDepth / mSampleMeanReadDepth;

        DoubleColumn gcNormalisedRatio = ratiosWithMedianCount.doubleColumn(CobaltColumns.RATIO)
                .multiply(medianNormalisation)
                .divide(ratiosWithMedianCount.doubleColumn("gcMedianCount"))
                .map(d -> Double.isFinite(d) ? d : Double.NaN); // protect against division by 0

        ratiosWithMedianCount.replaceColumn(gcNormalisedRatio.setName(CobaltColumns.RATIO));

        // resort it, the join messes up with the ordering
        ratiosWithMedianCount = ratiosWithMedianCount.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        mGCMedianReadDepth = gcMedianCalcDf;
        return ratiosWithMedianCount;
    }

    public Table gcMedianReadDepthTable()
    {
        return mGCMedianReadDepth;
    }

    public double getSampleMedianReadDepth()
    {
        return mSampleMedianReadDepth;
    }

    public double getSampleMeanReadDepth()
    {
        return mSampleMeanReadDepth;
    }

    // convert the gc median read count table to the object representation
    public GCMedianReadDepth gcMedianReadDepth()
    {
        final Map<GCBucket, Double> medianPerBucket = new HashMap<>();
        for (Row row : mGCMedianReadDepth)
        {
            medianPerBucket.put(new ImmutableGCBucket(row.getInt(CobaltColumns.GC_BUCKET)), row.getDouble("gcMedianCount"));
        }
        return new GCMedianReadDepth(mSampleMeanReadDepth, mSampleMedianReadDepth, medianPerBucket);
    }
}
