package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

// run-level counters for SpliceLiftBack: per-category, alignment-set composition × MAPQ tier, category × MAPQ tier.
// emits a summary at end-of-run (logs + a TSV) so reviewers can see prevalence before any category pruning is committed.
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

    public void record(final SAMRecord record, final LiftBackResult result)
    {
        // pre-drop full-set composition: includes alts dropped by the discriminator, so the summary reflects what
        // bwa actually produced (vs. result.comp(), the post-drop "kept" view written to TSV-A).
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

    public void logSummary()
    {
        RD_LOGGER.info("processed {} records", mTotal);

        RD_LOGGER.info("alignment-set composition x MAPQ tier:");
        logTable(mCompositionByTier, LiftBackResult.Composition.values(), MapqTier.values());

        RD_LOGGER.info("category x MAPQ tier:");
        logTable(mCategoryByTier, LiftBackCategory.values(), MapqTier.values());

        RD_LOGGER.info("primary bucket -> category breakdown:");
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
            RD_LOGGER.info("  {}: total={} {}", bucket.name(), bucketTotal, sb.toString());
        }

        int unliftable = mPerCategory[LiftBackCategory.LIFT_FAILED.ordinal()];
        if(mTotal > 0 && unliftable / (double) mTotal > LIFT_FAILED_WARN_FRACTION)
        {
            // run is already done; ERROR (not throw) so the operator sees this in log scrapers but the
            // BAM/TSVs we've already written are still available for inspection.
            RD_LOGGER.error("LIFT_FAILED rate {} / {} = {}% exceeds {}% threshold — likely sidecar/FASTA mismatch",
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

            // bucket -> category breakdown: retains the secondary (full) category alongside the
            // primary bucket grouping used for headline scanning.
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

        RD_LOGGER.info("wrote summary to {}", path);
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

            RD_LOGGER.info("  {}: total={} {}", row.name(), rowTotal, sb.toString());
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
