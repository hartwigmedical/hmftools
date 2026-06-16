package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

// Run-level counters for SpliceLiftBack: per-category, composition x MAPQ tier, category x MAPQ tier.
public class LiftBackStats
{
    private static final double LIFT_FAILED_WARN_FRACTION = 0.01;

    public enum MapqTier
    {
        MAPQ_ZERO,         // bwa flagged as multi-mapping
        MAPQ_POS_MULTI,    // MAPQ > 0 but XA non-empty (some ambiguity remains)
        MAPQ_POS_UNIQUE    // MAPQ > 0, no XA alts
    }

    private static final int N_CATEGORIES = LiftBackCategory.values().length;
    private static final int N_COMPOSITIONS = LiftBackResult.Composition.values().length;
    private static final int N_TIERS = MapqTier.values().length;

    private final int[] mPerCategory = new int[N_CATEGORIES];
    private final int[][] mCompositionByTier = new int[N_COMPOSITIONS][N_TIERS];
    private final int[][] mCategoryByTier = new int[N_CATEGORIES][N_TIERS];
    private int mTotal = 0;
    private int mLowAsSuppsDropped = 0;
    private int mOrphanSuppsDropped = 0;
    private int mLowAsPrimariesUnmapped = 0;

    public void record(final SAMRecord record, final LiftBackResult result)
    {
        // Full pre-drop composition (includes discriminator-dropped alts) reflects what bwa produced.
        LiftBackResult.Composition fullComposition = LiftBackResult.Composition.fromAlignments(result.liftedAlignments());
        MapqTier tier = deriveMapqTier(record, result);

        ++mTotal;
        ++mPerCategory[result.category().ordinal()];
        ++mCompositionByTier[fullComposition.ordinal()][tier.ordinal()];
        ++mCategoryByTier[result.category().ordinal()][tier.ordinal()];
    }

    public int total()
    {
        return mTotal;
    }

    public int categoryCount(final LiftBackCategory category)
    {
        return mPerCategory[category.ordinal()];
    }

    public void recordLowAsSuppDropped()
    {
        ++mLowAsSuppsDropped;
    }

    public int lowAsSuppsDropped()
    {
        return mLowAsSuppsDropped;
    }

    public void recordOrphanSuppDropped()
    {
        ++mOrphanSuppsDropped;
    }

    public int orphanSuppsDropped()
    {
        return mOrphanSuppsDropped;
    }

    public void recordLowAsPrimaryUnmapped()
    {
        ++mLowAsPrimariesUnmapped;
    }

    public int lowAsPrimariesUnmapped()
    {
        return mLowAsPrimariesUnmapped;
    }

    // Merge per-worker stats into this instance for a combined end-of-run summary.
    public void merge(final LiftBackStats other)
    {
        mTotal += other.mTotal;
        mLowAsSuppsDropped += other.mLowAsSuppsDropped;
        mOrphanSuppsDropped += other.mOrphanSuppsDropped;
        mLowAsPrimariesUnmapped += other.mLowAsPrimariesUnmapped;
        for(int i = 0; i < mPerCategory.length; ++i)
            mPerCategory[i] += other.mPerCategory[i];
        for(int i = 0; i < mCompositionByTier.length; ++i)
            for(int j = 0; j < N_TIERS; ++j)
                mCompositionByTier[i][j] += other.mCompositionByTier[i][j];
        for(int i = 0; i < mCategoryByTier.length; ++i)
            for(int j = 0; j < N_TIERS; ++j)
                mCategoryByTier[i][j] += other.mCategoryByTier[i][j];
    }

    public void logSummary()
    {
        TARS_LOGGER.info("processed {} records", mTotal);
        if(mLowAsSuppsDropped > 0)
            TARS_LOGGER.info("dropped {} non-rescued supps with source AS < {}",
                    mLowAsSuppsDropped, LiftBackGroupProcessor.SUPP_AS_DROP_THRESHOLD);
        if(mOrphanSuppsDropped > 0)
            TARS_LOGGER.info("dropped {} supps whose SA entries all failed to lift (orphaned in genomic space)",
                    mOrphanSuppsDropped);
        if(mLowAsPrimariesUnmapped > 0)
            TARS_LOGGER.info("unmapped {} primaries with AS < {} not improved by liftback",
                    mLowAsPrimariesUnmapped, LiftBackGroupProcessor.PRIMARY_AS_UNMAP_THRESHOLD);

        TARS_LOGGER.info("alignment-set composition x MAPQ tier:");
        logTable(mCompositionByTier, LiftBackResult.Composition.values(), MapqTier.values());

        TARS_LOGGER.info("category x MAPQ tier:");
        logTable(mCategoryByTier, LiftBackCategory.values(), MapqTier.values());

        TARS_LOGGER.info("primary bucket -> category breakdown:");
        for(LiftBackCategory.PrimaryBucket bucket : LiftBackCategory.PrimaryBucket.values())
        {
            int bucketTotal = 0;
            StringBuilder sb = new StringBuilder();
            for(LiftBackCategory cat : LiftBackCategory.values())
            {
                if(cat.primaryBucket() != bucket)
                    continue;
                int count = mPerCategory[cat.ordinal()];
                bucketTotal += count;
                if(count > 0)
                    sb.append(cat.name()).append("=").append(count).append(" ");
            }
            if(bucketTotal == 0)
                continue;
            TARS_LOGGER.info("  {}: total={} {}", bucket.name(), bucketTotal, sb.toString());
        }

        int unliftable = mPerCategory[LiftBackCategory.LIFT_FAILED.ordinal()];
        if(mTotal > 0 && unliftable / (double) mTotal > LIFT_FAILED_WARN_FRACTION)
        {
            // ERROR not throw: BAM/TSVs already written should remain available for inspection.
            TARS_LOGGER.error("LIFT_FAILED rate {} / {} = {}% exceeds {}% threshold - likely sidecar/FASTA mismatch",
                    unliftable, mTotal,
                    String.format("%.2f", 100.0 * unliftable / mTotal),
                    String.format("%.2f", 100.0 * LIFT_FAILED_WARN_FRACTION));
        }
    }

    public void writeSummary(final String path) throws IOException
    {
        try(BufferedWriter writer = createBufferedWriter(path))
        {
            writer.write(String.join(TSV_DELIM, "Section", "RowKey", "ColumnKey", "Count"));
            writer.newLine();

            for(LiftBackResult.Composition comp : LiftBackResult.Composition.values())
            {
                for(MapqTier tier : MapqTier.values())
                {
                    int count = mCompositionByTier[comp.ordinal()][tier.ordinal()];
                    if(count == 0)
                        continue;
                    writer.write(String.join(TSV_DELIM,
                            "composition_x_mapq", comp.name(), tier.name(), String.valueOf(count)));
                    writer.newLine();
                }
            }

            for(LiftBackCategory cat : LiftBackCategory.values())
            {
                for(MapqTier tier : MapqTier.values())
                {
                    int count = mCategoryByTier[cat.ordinal()][tier.ordinal()];
                    if(count == 0)
                        continue;
                    writer.write(String.join(TSV_DELIM,
                            "category_x_mapq", cat.name(), tier.name(), String.valueOf(count)));
                    writer.newLine();
                }
            }

            for(LiftBackCategory.PrimaryBucket bucket : LiftBackCategory.PrimaryBucket.values())
            {
                for(LiftBackCategory cat : LiftBackCategory.values())
                {
                    if(cat.primaryBucket() != bucket)
                        continue;
                    int count = mPerCategory[cat.ordinal()];
                    if(count == 0)
                        continue;
                    writer.write(String.join(TSV_DELIM,
                            "bucket_x_category", bucket.name(), cat.name(), String.valueOf(count)));
                    writer.newLine();
                }
            }
        }

        TARS_LOGGER.info("wrote summary to {}", path);
    }

    private static <R extends Enum<R>, C extends Enum<C>> void logTable(
            final int[][] table, final R[] rows, final C[] cols)
    {
        for(R row : rows)
        {
            StringBuilder sb = new StringBuilder();
            int rowTotal = 0;
            for(C col : cols)
            {
                int count = table[row.ordinal()][col.ordinal()];
                rowTotal += count;
                if(count > 0)
                    sb.append(col.name()).append("=").append(count).append(" ");
            }

            if(rowTotal == 0)
                continue;

            TARS_LOGGER.info("  {}: total={} {}", row.name(), rowTotal, sb.toString());
        }
    }

    static MapqTier deriveMapqTier(final SAMRecord record, final LiftBackResult result)
    {
        int mapq = record.getMappingQuality();
        if(mapq == 0)
            return MapqTier.MAPQ_ZERO;
        return result.numXaAlts() > 0 ? MapqTier.MAPQ_POS_MULTI : MapqTier.MAPQ_POS_UNIQUE;
    }
}
