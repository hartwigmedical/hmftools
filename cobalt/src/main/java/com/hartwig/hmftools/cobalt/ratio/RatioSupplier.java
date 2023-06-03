package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltUtils.toCommonChromosomeMap;
import static com.hartwig.hmftools.cobalt.ratio.DiploidRatioSupplier.calcDiploidRatioResults;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import tech.tablesaw.api.*;
import tech.tablesaw.columns.Column;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.diploid.DiploidRatioLoader;
import com.hartwig.hmftools.cobalt.lowcov.LowCovBucket;
import com.hartwig.hmftools.cobalt.lowcov.LowCoverageRatioBuilder;
import com.hartwig.hmftools.cobalt.targeted.TargetedRatioBuilder;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.MedianRatioFile;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCountFile;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class RatioSupplier
{
    private final String mTumorId;
    private final String mReferenceId;
    private final String mOutputDir;

    private final Collection<Chromosome> mChromosomes;
    private final Table mGcProfiles;
    @Nullable private final Table mReferenceCounts;
    @Nullable private final Table mTumorCounts;

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
        protected final GcNormalizedRatioBuilder gcNormalizedRatioBuilder;

        @Nullable Multimap<String, LowCovBucket> consolidatedBuckets;

        // chromosomePositionIndex, ratio
        Table readRatios;

        Table getRatios() { return readRatios; }

        SampleRatios(
                final String sampleId,
                final Table readCounts,
                final Table gcProfiles,
                @Nullable Table targetRegionEnrichment,
                SparseBucketPolicy sparseBucketPolicy,
                @Nullable Multimap<String, LowCovBucket> consolidatedBuckets,
                final String outputDir,
                ChromosomePositionCodec chromosomePosCodec) throws IOException
        {
            CB_LOGGER.info("calculating sample ratios for {}", sampleId);

            // just convert read count to ratios

            // We start by setting ratio as the read count
            readRatios = readCounts.copy();
            readRatios.addColumns(readRatios.intColumn(CobaltColumns.READ_COUNT).asDoubleColumn().setName(CobaltColumns.RATIO));

            CB_LOGGER.info("merging in gc profile");

            // merge in the gc profile
            readRatios = readRatios.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).leftOuter(gcProfiles);

            // on target ratios
            if (targetRegionEnrichment != null)
            {
                CB_LOGGER.info("using targeted ratio");
                readRatios = new TargetedRatioBuilder(readRatios, targetRegionEnrichment, chromosomePosCodec).ratios();
            }

            gcNormalizedRatioBuilder = new GcNormalizedRatioBuilder(readRatios, true);
            readRatios = gcNormalizedRatioBuilder.ratios();

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
                    double medianReadCount = gcNormalizedRatioBuilder.getSampleMedianReadCount();
                    this.consolidatedBuckets = LowCoverageRatioBuilder.calcConsolidateBuckets(readRatios, medianReadCount);
                    break;
                }
            }

            if (this.consolidatedBuckets != null)
            {
                CB_LOGGER.info("using low coverage ratio");
                readRatios = new LowCoverageRatioBuilder(readRatios, this.consolidatedBuckets, chromosomePosCodec).ratios();
            }

            CB_LOGGER.info("Persisting {} gc read count to {}", sampleId, outputDir);
            final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDir, sampleId);
            GCMedianReadCountFile.write(tumorGCMedianFilename, gcNormalizedRatioBuilder.gcMedianReadCount());
        }
    }

    static class GermlineRatios extends SampleRatios
    {
        // processing states
        private final Table gcDiploidRatios;

        GermlineRatios(final String referenceId,
                final Table readCounts,
                final Table gcProfiles,
                @Nullable Table targetRegionEnrichment,
                SparseBucketPolicy sparseBucketPolicy,
                @Nullable Multimap<String, LowCovBucket> consolidatedBuckets,
                final Collection<Chromosome> chromosomes,
                final String outputDir,
                ChromosomePositionCodec chromosomePosCodec) throws IOException
        {
            super(referenceId, readCounts, gcProfiles, targetRegionEnrichment, sparseBucketPolicy,
                    consolidatedBuckets, outputDir, chromosomePosCodec);

            // TODO: check this
            final List<MedianRatio> medianRatios = MedianRatioFactory.createFromReadRatio(toCommonChromosomeMap(getRatios()));

            CB_LOGGER.info("Persisting {} gc ratio medians to {}", referenceId, outputDir);
            final String ratioMedianFilename = MedianRatioFile.generateFilename(outputDir, referenceId);
            MedianRatioFile.write(ratioMedianFilename, medianRatios);

            CB_LOGGER.info("Applying ratio diploid normalization");
            gcDiploidRatios = calcDiploidRatioResults(chromosomes, getRatios(), medianRatios);
        }
    }

    public RatioSupplier(final String reference, final String tumor, final String outputDirectory,
            final Table gcProfiles,
            final Collection<Chromosome> chromosomes,
            @Nullable final Table referenceCounts,
            @Nullable final Table tumorCounts,
            ChromosomePositionCodec chromosomePosCodec)
    {
        mTumorId = tumor;
        mReferenceId = reference;
        mOutputDir = outputDirectory;
        mGcProfiles = gcProfiles;
        mChromosomes = chromosomes;
        mReferenceCounts = referenceCounts;
        mTumorCounts = tumorCounts;
        mChromosomePosCodec = chromosomePosCodec;
    }
    
    public void setTargetRegionEnrichment(Table targetRegionEnrichment)
    {
        mTargetRegionEnrichment = targetRegionEnrichment;
    }

    @NotNull
    public Table tumorOnly(final String diploidBedFile) throws IOException
    {
        if (mTumorCounts == null)
        {
            CB_LOGGER.fatal("Tumor count should not be null");
            throw new RuntimeException("tumor count is null");
        }
        SparseBucketPolicy sparseBucketPolicy = mTargetRegionEnrichment == null ? SparseBucketPolicy.CALC_CONSOLIDATED_BUCKETS : SparseBucketPolicy.DO_NOT_CONSOLIDATE;
        var tumorRatios = new SampleRatios(mTumorId, mTumorCounts, mGcProfiles, mTargetRegionEnrichment, sparseBucketPolicy,
                null, mOutputDir, mChromosomePosCodec);
        final Table diploidRatios = new DiploidRatioLoader(mChromosomes, diploidBedFile, mChromosomePosCodec).build();

        // merge this ratios together into one cobalt ratio
        return mergeRatios(
                null, mTumorCounts,
                diploidRatios, tumorRatios.getRatios(), diploidRatios);
    }

    @NotNull
    public Table germlineOnly() throws IOException
    {
        if (mReferenceCounts == null)
        {
            CB_LOGGER.fatal("Reference count should not be null");
            throw new RuntimeException("reference count is null");
        }
        SparseBucketPolicy sparseBucketPolicy = mTargetRegionEnrichment == null ? SparseBucketPolicy.CALC_CONSOLIDATED_BUCKETS : SparseBucketPolicy.DO_NOT_CONSOLIDATE;
        var germlineRatios = new GermlineRatios(mReferenceId, mReferenceCounts, mGcProfiles, mTargetRegionEnrichment,
                sparseBucketPolicy, null, mChromosomes, mOutputDir, mChromosomePosCodec);
        return mergeRatios(
                mReferenceCounts, null,
                germlineRatios.getRatios(), null, germlineRatios.gcDiploidRatios);
    }

    @NotNull
    public Table tumorNormalPair() throws IOException
    {
        if (mReferenceCounts == null)
        {
            CB_LOGGER.fatal("Reference count should not be null");
            throw new RuntimeException("reference count is null");
        }
        if (mTumorCounts == null)
        {
            CB_LOGGER.fatal("Tumor count should not be null");
            throw new RuntimeException("tumor count is null");
        }
        SparseBucketPolicy tumorSparseBucketPolicy = mTargetRegionEnrichment == null ?
                SparseBucketPolicy.CALC_CONSOLIDATED_BUCKETS : SparseBucketPolicy.DO_NOT_CONSOLIDATE;

        var tumorRatios = new SampleRatios(mTumorId, mTumorCounts, mGcProfiles, mTargetRegionEnrichment,
                tumorSparseBucketPolicy, null, mOutputDir, mChromosomePosCodec);

        SparseBucketPolicy germlineSparseBucketPolicy = tumorRatios.consolidatedBuckets == null ?
                SparseBucketPolicy.DO_NOT_CONSOLIDATE : SparseBucketPolicy.USE_PROVIDED_BUCKETS;

        var germlineRatios = new GermlineRatios(mReferenceId, mReferenceCounts, mGcProfiles, mTargetRegionEnrichment,
                germlineSparseBucketPolicy, tumorRatios.consolidatedBuckets, mChromosomes, mOutputDir, mChromosomePosCodec);

        return mergeRatios(
                mReferenceCounts, mTumorCounts,
                germlineRatios.getRatios(), tumorRatios.getRatios(), germlineRatios.gcDiploidRatios);
    }

    // merge everything together
    @NotNull
    private static Table mergeRatios(
            @Nullable Table referenceCounts,
            @Nullable Table tumorCounts,
            @Nullable Table referenceRatios,
            @Nullable Table tumorRatios,
            @Nullable Table referenceDiploidRatios)
    {
        CB_LOGGER.info("merging ratios");

        // get all the chromosome positions from the counts
        Table result = Table.create(LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS));

        // now we make sure all tables are valid
        if (referenceCounts == null)
        {
            referenceCounts = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    IntColumn.create(CobaltColumns.READ_COUNT));
        }

        if (tumorCounts == null)
        {
            tumorCounts = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    IntColumn.create(CobaltColumns.READ_COUNT));
        }

        if (referenceRatios == null)
        {
            referenceRatios = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    IntColumn.create(CobaltColumns.RATIO));
        }

        if (tumorRatios == null)
        {
            tumorRatios = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    IntColumn.create(CobaltColumns.RATIO));
        }

        if (referenceDiploidRatios == null)
        {
            referenceDiploidRatios = Table.create(
                    LongColumn.create(CobaltColumns.ENCODED_CHROMOSOME_POS),
                    IntColumn.create(CobaltColumns.RATIO));
        }

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                referenceCounts.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.READ_COUNT));

        // rename the readCount column
        result.intColumn("readCount").setName("referenceReadCount");

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                tumorCounts.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.READ_COUNT));

        // rename the readCount column
        result.intColumn("readCount").setName("tumorReadCount");

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                referenceRatios.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.RATIO));

        // rename the ratio column
        result.doubleColumn("ratio").setName("referenceGCRatio");

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                tumorRatios.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.RATIO));

        // rename the ratio column
        result.doubleColumn("ratio").setName("tumorGCRatio");

        result = result.joinOn(CobaltColumns.ENCODED_CHROMOSOME_POS).fullOuter(
                referenceDiploidRatios.retainColumns(CobaltColumns.ENCODED_CHROMOSOME_POS, CobaltColumns.RATIO));

        // rename the ratio column
        result.doubleColumn("ratio").setName("referenceGCDiploidRatio");

        Validate.isTrue(result.longColumn(CobaltColumns.ENCODED_CHROMOSOME_POS).isMissing().isEmpty());

        // sort by the encoded chromosome pos
        result = result.sortAscendingOn(CobaltColumns.ENCODED_CHROMOSOME_POS);

        // set any missing value to -1
        for (Column<?> c: result.columns())
        {
            if (c instanceof IntColumn)
            {
                ((IntColumn)c).setMissingTo(-1);
            }
            else if (c instanceof DoubleColumn)
            {
                ((DoubleColumn)c).setMissingTo(-1.0);
            }
        }

        return result;
    }
}