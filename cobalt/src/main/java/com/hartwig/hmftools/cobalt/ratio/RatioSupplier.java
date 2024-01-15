package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltUtils.toCommonChromosomeMap;
import static com.hartwig.hmftools.cobalt.ratio.DiploidRatioSupplier.calcDiploidRatioResults;

import java.io.IOException;
import java.util.List;

import tech.tablesaw.api.*;
import tech.tablesaw.columns.Column;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.cobalt.lowcov.LowCovBucket;
import com.hartwig.hmftools.cobalt.lowcov.LowCoverageRatioMapper;
import com.hartwig.hmftools.cobalt.targeted.TargetedRatioMapper;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatioFile;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadDepthFile;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class RatioSupplier
{
    private final String mTumorId;
    private final String mReferenceId;
    @Nullable private final String mOutputDir;

    private final Table mGcProfiles;
    @Nullable private final Table mReferenceDepths;
    @Nullable private final Table mTumorDepths;

    // a table with chromosome, position, relativeEnrichment
    private Table mTargetRegionEnrichment = null;

    private final ChromosomePositionCodec mChromosomePosCodec;

    enum SparseBucketPolicy
    {
        DO_NOT_CONSOLIDATE,
        USE_PROVIDED_BUCKETS,
        CALC_CONSOLIDATED_BUCKETS
    }

    static class SampleRatios
    {
        // processing states
        protected final GcNormalizedRatioMapper gcNormalizedRatioMapper;

        @Nullable Multimap<String, LowCovBucket> consolidatedBuckets;

        // chromosomePositionIndex, ratio
        Table readRatios;

        Table getRatios() { return readRatios; }

        SampleRatios(
                final String sampleId,
                final Table readDepths,
                final Table gcProfiles,
                @Nullable Table targetRegionEnrichment,
                SparseBucketPolicy sparseBucketPolicy,
                @Nullable Multimap<String, LowCovBucket> consolidatedBuckets,
                @Nullable final String outputDir,
                ChromosomePositionCodec chromosomePosCodec) throws IOException
        {
            CB_LOGGER.info("calculating sample ratios for {}", sampleId);

            readRatios = readDepths.copy().setName("readRatios");

            CB_LOGGER.info("merging in GC profile");

            // merge in the gc profile
            readRatios = readRatios.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).leftOuter(gcProfiles);

            // set column as ratio, but filter out unmappable regions
            DoubleColumn ratioColumn = DoubleColumn.create(CobaltColumns.RATIO);
            for (Row row : readRatios)
            {
                if(!row.isMissing(CobaltColumns.IS_MAPPABLE) && row.getBoolean(CobaltColumns.IS_MAPPABLE))
                {
                    ratioColumn.append(row.getDouble(CobaltColumns.READ_DEPTH));
                }
                else
                {
                    ratioColumn.appendMissing();
                }
            }
            readRatios.addColumns(ratioColumn);

            // on target ratios
            if(targetRegionEnrichment != null)
            {
                CB_LOGGER.info("using targeted ratio");
                readRatios = new TargetedRatioMapper(targetRegionEnrichment, chromosomePosCodec).mapRatios(readRatios);
            }

            gcNormalizedRatioMapper = new GcNormalizedRatioMapper();
            readRatios = gcNormalizedRatioMapper.mapRatios(readRatios);

            switch (sparseBucketPolicy)
            {
                case DO_NOT_CONSOLIDATE:
                    this.consolidatedBuckets = null;
                    break;
                case USE_PROVIDED_BUCKETS:
                    this.consolidatedBuckets = consolidatedBuckets;
                    break;
                case CALC_CONSOLIDATED_BUCKETS:
                {
                    // determine consolidated buckets
                    // determine the low cov consolidation window count
                    double medianReadDepth = gcNormalizedRatioMapper.getSampleMedianReadDepth();
                    this.consolidatedBuckets = LowCoverageRatioMapper.calcConsolidateBuckets(readRatios, medianReadDepth);
                    break;
                }
            }

            if(this.consolidatedBuckets != null)
            {
                CB_LOGGER.info("using low coverage ratio");
                readRatios = new LowCoverageRatioMapper(this.consolidatedBuckets, chromosomePosCodec).mapRatios(readRatios);
            }

            if(outputDir != null)
            {
                CB_LOGGER.info("persisting {} gc read count to {}", sampleId, outputDir);
                final String tumorGCMedianFilename = GCMedianReadDepthFile.generateFilename(outputDir, sampleId);
                GCMedianReadDepthFile.write(tumorGCMedianFilename, gcNormalizedRatioMapper.gcMedianReadDepth());
            }
        }
    }

    static class GermlineRatios extends SampleRatios
    {
        // processing states
        private final Table gcDiploidRatios;

        GermlineRatios(final String referenceId,
                final Table readDepths,
                final Table gcProfiles,
                @Nullable Table targetRegionEnrichment,
                SparseBucketPolicy sparseBucketPolicy,
                @Nullable Multimap<String, LowCovBucket> consolidatedBuckets,
                final String outputDir,
                ChromosomePositionCodec chromosomePosCodec) throws IOException
        {
            super(referenceId, readDepths, gcProfiles, targetRegionEnrichment, sparseBucketPolicy,
                    consolidatedBuckets, outputDir, chromosomePosCodec);

            // TODO: check this
            final List<MedianRatio> medianRatios = MedianRatioFactory.createFromReadRatio(toCommonChromosomeMap(getRatios()));

            CB_LOGGER.info("persisting {} gc ratio medians to {}", referenceId, outputDir);
            final String ratioMedianFilename = MedianRatioFile.generateFilename(outputDir, referenceId);
            MedianRatioFile.write(ratioMedianFilename, medianRatios);

            CB_LOGGER.info("applying ratio diploid normalization");
            gcDiploidRatios = calcDiploidRatioResults(getRatios(), medianRatios);
        }
    }

    public RatioSupplier(final String reference, final String tumor,
            @Nullable final String outputDirectory,
            final Table gcProfiles,
            @Nullable final Table referenceDepths,
            @Nullable final Table tumorDepths,
            ChromosomePositionCodec chromosomePosCodec)
    {
        mTumorId = tumor;
        mReferenceId = reference;
        mOutputDir = outputDirectory;
        mGcProfiles = gcProfiles;
        mReferenceDepths = referenceDepths;
        mTumorDepths = tumorDepths;
        mChromosomePosCodec = chromosomePosCodec;
    }
    
    public void setTargetRegionEnrichment(Table targetRegionEnrichment)
    {
        mTargetRegionEnrichment = targetRegionEnrichment;
    }

    @NotNull
    public Table tumorOnly(final Table diploidRegions) throws IOException
    {
        if(mTumorDepths == null)
        {
            CB_LOGGER.error("tumor count should not be null");
            throw new RuntimeException("tumor count is null");
        }
        SparseBucketPolicy sparseBucketPolicy = mTargetRegionEnrichment == null ? SparseBucketPolicy.CALC_CONSOLIDATED_BUCKETS : SparseBucketPolicy.DO_NOT_CONSOLIDATE;
        Table tumorRatios = new SampleRatios(mTumorId, mTumorDepths, mGcProfiles, mTargetRegionEnrichment, sparseBucketPolicy,
                null, mOutputDir, mChromosomePosCodec).getRatios();

        // filter tumor ratios by the diploid regions
        // we use inner join to remove any tumor ratios that are not in the diploid regions
        tumorRatios = tumorRatios.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS)
                .inner(diploidRegions.selectColumns(CobaltColumns.ENCODED_CHROMOSOME_POS))
                .sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        // merge this ratios together into one cobalt ratio
        return mergeRatios(null, mTumorDepths, null, tumorRatios, null);
    }

    @NotNull
    public Table germlineOnly() throws IOException
    {
        if(mReferenceDepths == null)
        {
            CB_LOGGER.fatal("Reference count should not be null");
            throw new RuntimeException("reference count is null");
        }
        SparseBucketPolicy sparseBucketPolicy = mTargetRegionEnrichment == null ? SparseBucketPolicy.CALC_CONSOLIDATED_BUCKETS : SparseBucketPolicy.DO_NOT_CONSOLIDATE;
        var germlineRatios = new GermlineRatios(mReferenceId, mReferenceDepths, mGcProfiles, mTargetRegionEnrichment,
                sparseBucketPolicy, null, mOutputDir, mChromosomePosCodec);
        return mergeRatios(
                mReferenceDepths, null,
                germlineRatios.getRatios(), null, germlineRatios.gcDiploidRatios);
    }

    @NotNull
    public Table tumorNormalPair() throws IOException
    {
        if(mReferenceDepths == null)
        {
            CB_LOGGER.fatal("Reference count should not be null");
            throw new RuntimeException("reference count is null");
        }
        if(mTumorDepths == null)
        {
            CB_LOGGER.fatal("Tumor count should not be null");
            throw new RuntimeException("tumor count is null");
        }
        SparseBucketPolicy tumorSparseBucketPolicy = mTargetRegionEnrichment == null ?
                SparseBucketPolicy.CALC_CONSOLIDATED_BUCKETS : SparseBucketPolicy.DO_NOT_CONSOLIDATE;

        var tumorRatios = new SampleRatios(mTumorId, mTumorDepths, mGcProfiles, mTargetRegionEnrichment,
                tumorSparseBucketPolicy, null, mOutputDir, mChromosomePosCodec);

        SparseBucketPolicy germlineSparseBucketPolicy = tumorRatios.consolidatedBuckets == null ?
                SparseBucketPolicy.DO_NOT_CONSOLIDATE : SparseBucketPolicy.USE_PROVIDED_BUCKETS;

        var germlineRatios = new GermlineRatios(mReferenceId, mReferenceDepths, mGcProfiles, mTargetRegionEnrichment,
                germlineSparseBucketPolicy, tumorRatios.consolidatedBuckets, mOutputDir, mChromosomePosCodec);

        return mergeRatios(
                mReferenceDepths, mTumorDepths,
                germlineRatios.getRatios(), tumorRatios.getRatios(), germlineRatios.gcDiploidRatios);
    }

    // merge everything together
    @NotNull
    private static Table mergeRatios(
            @Nullable Table referenceDepths,
            @Nullable Table tumorDepths,
            @Nullable Table referenceRatios,
            @Nullable Table tumorRatios,
            @Nullable Table referenceDiploidRatios)
    {
        CB_LOGGER.info("start merging ratios");

        // get all the chromosome positions from the counts
        Table result = Table.create(LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS));

        // now we make sure all tables are valid, by setting missing tables to empty
        if(referenceDepths == null)
        {
            referenceDepths = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    DoubleColumn.create(CobaltColumns.READ_DEPTH),
                    DoubleColumn.create(CobaltColumns.READ_GC_CONTENT));
        }

        if(tumorDepths == null)
        {
            tumorDepths = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    DoubleColumn.create(CobaltColumns.READ_DEPTH),
                    DoubleColumn.create(CobaltColumns.READ_GC_CONTENT));
        }

        if(referenceRatios == null)
        {
            referenceRatios = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    DoubleColumn.create(CobaltColumns.RATIO));
        }

        if(tumorRatios == null)
        {
            tumorRatios = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    DoubleColumn.create(CobaltColumns.RATIO));
        }

        if(referenceDiploidRatios == null)
        {
            referenceDiploidRatios = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    DoubleColumn.create(CobaltColumns.RATIO));
        }

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                referenceDepths.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.READ_DEPTH, CobaltColumns.READ_GC_CONTENT));

        // rename the readDepth column
        result.doubleColumn(CobaltColumns.READ_DEPTH).setName(CobaltColumns.REFERENCE_READ_DEPTH);
        result.doubleColumn(CobaltColumns.READ_GC_CONTENT).setName(CobaltColumns.REFERENCE_GC_CONTENT);

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                tumorDepths.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.READ_DEPTH, CobaltColumns.READ_GC_CONTENT));

        // rename the readDepth column
        result.doubleColumn(CobaltColumns.READ_DEPTH).setName(CobaltColumns.TUMOR_READ_DEPTH);
        result.doubleColumn(CobaltColumns.READ_GC_CONTENT).setName(CobaltColumns.TUMOR_GC_CONTENT);

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                referenceRatios.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.RATIO));

        // rename the ratio column
        result.doubleColumn(CobaltColumns.RATIO).setName("referenceGCRatio");

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                tumorRatios.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.RATIO));

        // rename the ratio column
        result.doubleColumn(CobaltColumns.RATIO).setName("tumorGCRatio");

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                referenceDiploidRatios.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.RATIO));

        // rename the ratio column
        result.doubleColumn(CobaltColumns.RATIO).setName("referenceGCDiploidRatio");

        Validate.isTrue(result.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).isMissing().isEmpty());

        // sort by the encoded chromosome pos
        result = result.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        // set any missing value to -1
        for (Column<?> c: result.columns())
        {
            if(c instanceof IntColumn)
            {
                ((IntColumn)c).setMissingTo(CobaltConstants.INVALID_VALUE_INDICATOR);
            }
            else if(c instanceof DoubleColumn)
            {
                ((DoubleColumn)c).setMissingTo((double)CobaltConstants.INVALID_VALUE_INDICATOR);
                // ((DoubleColumn)c).map(d -> Double.isFinite(d) ? d : -1.0);
            }
        }

        CB_LOGGER.info("finish merging ratios");

        return result;
    }
}