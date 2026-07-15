package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_AS_DROP_THRESHOLD;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

// Run-level counters for TarsApplication: per deciding-feature, per record-state, composition x MAPQ tier,
// deciding-feature x MAPQ tier. Non-RESOLVED records (unmapped / lift-failed / supplementary) carry no
// deciding feature, so they are counted only in the record-state tally.
public class LiftBackStats
{
    private static final double LIFT_FAILED_WARN_FRACTION = 0.01;

    public enum MapqTier
    {
        MAPQ_ZERO,         // bwa flagged as multi-mapping
        MAPQ_POS_MULTI,    // MAPQ > 0 but XA non-empty (some ambiguity remains)
        MAPQ_POS_UNIQUE    // MAPQ > 0, no XA alts
    }

    private static final int N_FEATURES = DecidingFeature.values().length;
    private static final int N_STATES = RecordState.values().length;
    private static final int N_COMPOSITIONS = LiftBackResult.Composition.values().length;
    private static final int N_TIERS = MapqTier.values().length;

    private final int[] mPerFeature = new int[N_FEATURES];
    private final int[] mSwapByFeature = new int[N_FEATURES];
    private final int[] mPerState = new int[N_STATES];
    private final int[][] mCompositionByTier = new int[N_COMPOSITIONS][N_TIERS];
    private final int[][] mFeatureByTier = new int[N_FEATURES][N_TIERS];
    private int mTotal = 0;
    private int mSplicedOutput = 0;        // resolved records emitted with an N (spliced)
    private int mMapqZeroIn = 0;           // bwa MAPQ 0 on a resolved primary
    private int mMapqRescued = 0;          // of those, lifted to a positive MAPQ
    private int mLowAsSuppsDropped = 0;
    private int mOrphanSuppsDropped = 0;
    private int mLowAsPrimariesUnmapped = 0;
    private PassEffects mPassEffects = PassEffects.EMPTY;

    // Per-pass effectiveness, aggregated across workers by the caller and folded into the summary so the
    // always-written summary answers "what did liftback do" without the opt-in ~100GB per-record TSV.
    public record PassEffects(
            int suppCandidates, int suppMerged, int suppClamped,
            long collapseLeading, long collapseTrailing,
            long reclaimRecords, long reclaimBasesLead, long reclaimBasesTrail, long altsDropped,
            long overCapUnmapped, long excludedReads)
    {
        public static final PassEffects EMPTY = new PassEffects(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }

    public void setPassEffects(final PassEffects passEffects)
    {
        mPassEffects = passEffects;
    }

    public void record(final SAMRecord record, final LiftBackResult result)
    {
        // Full pre-drop composition (includes discriminator-dropped alts) reflects what bwa produced.
        LiftBackResult.Composition fullComposition = LiftBackResult.Composition.fromAlignments(result.liftedAlignments());
        MapqTier tier = deriveMapqTier(record, result);

        ++mTotal;
        ++mCompositionByTier[fullComposition.ordinal()][tier.ordinal()];
        ++mPerState[result.recordState().ordinal()];

        DecidingFeature feature = result.decidingFeature();
        if(feature != null) // RESOLVED records only
        {
            ++mPerFeature[feature.ordinal()];
            ++mFeatureByTier[feature.ordinal()][tier.ordinal()];
            if(result.swapped())
            {
                ++mSwapByFeature[feature.ordinal()];
            }
            if(result.hasNCigar())
            {
                ++mSplicedOutput;
            }
            if(result.inputMapq() == 0)
            {
                ++mMapqZeroIn;
                if(result.updatedMapq() > 0)
                {
                    ++mMapqRescued;
                }
            }
        }
    }

    public int total()
    {
        return mTotal;
    }

    public int featureCount(final DecidingFeature feature)
    {
        return mPerFeature[feature.ordinal()];
    }

    public int stateCount(final RecordState state)
    {
        return mPerState[state.ordinal()];
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
        mSplicedOutput += other.mSplicedOutput;
        mMapqZeroIn += other.mMapqZeroIn;
        mMapqRescued += other.mMapqRescued;
        mLowAsSuppsDropped += other.mLowAsSuppsDropped;
        mOrphanSuppsDropped += other.mOrphanSuppsDropped;
        mLowAsPrimariesUnmapped += other.mLowAsPrimariesUnmapped;
        for(int i = 0; i < N_FEATURES; ++i)
        {
            mPerFeature[i] += other.mPerFeature[i];
            mSwapByFeature[i] += other.mSwapByFeature[i];
        }
        for(int i = 0; i < N_STATES; ++i)
        {
            mPerState[i] += other.mPerState[i];
        }
        for(int i = 0; i < mCompositionByTier.length; ++i)
            for(int j = 0; j < N_TIERS; ++j)
            {
                mCompositionByTier[i][j] += other.mCompositionByTier[i][j];
            }
        for(int i = 0; i < N_FEATURES; ++i)
            for(int j = 0; j < N_TIERS; ++j)
            {
                mFeatureByTier[i][j] += other.mFeatureByTier[i][j];
            }
    }

    // Outcome counts are derived from the per-feature tallies (outcome is a function of the deciding feature),
    // so they cover resolved primaries only - matching the old per-outcome counter.
    private int outcomeCount(final Outcome outcome)
    {
        int total = 0;
        for(DecidingFeature feature : DecidingFeature.values())
        {
            if(feature.outcome() == outcome)
            {
                total += mPerFeature[feature.ordinal()];
            }
        }
        return total;
    }

    public void logSummary()
    {
        TARS_LOGGER.info("processed {} records", mTotal);
        TARS_LOGGER.info("resolved outcomes: REF={} TX={} UNRESOLVED={}",
                outcomeCount(Outcome.REF), outcomeCount(Outcome.TX), outcomeCount(Outcome.UNRESOLVED));
        TARS_LOGGER.info("spliced output reads: {}; MAPQ-0 in: {}, rescued: {}", mSplicedOutput, mMapqZeroIn, mMapqRescued);
        if(mLowAsSuppsDropped > 0)
        {
            TARS_LOGGER.info("dropped {} non-resolved supps with source AS < {}",
                    mLowAsSuppsDropped, SUPP_AS_DROP_THRESHOLD);
        }
        if(mOrphanSuppsDropped > 0)
        {
            TARS_LOGGER.info("dropped {} supps whose SA entries all failed to lift (orphaned in genomic space)",
                    mOrphanSuppsDropped);
        }
        if(mLowAsPrimariesUnmapped > 0)
        {
            TARS_LOGGER.info("unmapped {} primaries with AS < {} not improved by liftback",
                    mLowAsPrimariesUnmapped, PRIMARY_AS_UNMAP_THRESHOLD);
        }

        TARS_LOGGER.debug("alignment-set composition x MAPQ tier:");
        logTable(mCompositionByTier, LiftBackResult.Composition.values(), MapqTier.values());

        TARS_LOGGER.debug("deciding-feature x MAPQ tier:");
        logTable(mFeatureByTier, DecidingFeature.values(), MapqTier.values());

        TARS_LOGGER.debug("record-state breakdown:");
        for(RecordState state : RecordState.values())
        {
            int count = mPerState[state.ordinal()];
            if(count > 0)
            {
                TARS_LOGGER.debug("  {}={}", state.name(), count);
            }
        }

        int unliftable = mPerState[RecordState.LIFT_FAILED.ordinal()];
        if(mTotal > 0 && unliftable / (double) mTotal > LIFT_FAILED_WARN_FRACTION)
        {
            // ERROR not throw: BAM/TSVs already written should remain available for inspection.
            TARS_LOGGER.error("LIFT_FAILED rate {} / {} = {}% exceeds {}% threshold - likely sidecar/FASTA mismatch",
                    unliftable, mTotal,
                    String.format("%.2f", 100.0 * unliftable / mTotal),
                    String.format("%.2f", 100.0 * LIFT_FAILED_WARN_FRACTION));
        }
    }

    // Flat, greppable one-metric-per-line summary: Metric, Value, Pct, Basis. Pct is Value as a percentage
    // of the named Basis metric so each rate is self-documenting (grep the metric, read the number and what
    // it is a fraction of). The 2D cross-tabs (composition/feature x MAPQ tier) stay in the DEBUG log only.
    public void writeSummary(final String path) throws IOException
    {
        long resolved = mPerState[RecordState.RESOLVED.ordinal()];
        long liftFailed = mPerState[RecordState.LIFT_FAILED.ordinal()];
        long swaps = 0;
        for(int i = 0; i < N_FEATURES; ++i)
        {
            swaps += mSwapByFeature[i];
        }

        try(BufferedWriter writer = createBufferedWriter(path))
        {
            writer.write(String.join(TSV_DELIM, "Metric", "Value", "Pct", "Basis"));
            writer.newLine();

            // records read, the primaries the discriminator placed (the denominator for the splits below), and
            // the lift-failure count - a non-zero rate here flags a sidecar/FASTA mismatch.
            writeMetric(writer, "records_total", mTotal);
            writeMetric(writer, "primaries_resolved", resolved, "records_total", mTotal);
            writeMetric(writer, "lift_failed", liftFailed, "records_total", mTotal);

            // discriminator outcome and the feature that decided it, over resolved primaries.
            for(Outcome outcome : Outcome.values())
            {
                writeMetric(writer, "outcome_" + outcome.name().toLowerCase(), outcomeCount(outcome), "primaries_resolved", resolved);
            }
            for(DecidingFeature feature : DecidingFeature.values())
            {
                writeMetric(
                        writer, "feature_" + feature.name().toLowerCase(), mPerFeature[feature.ordinal()],
                        "primaries_resolved", resolved);
            }
            writeMetric(writer, "primary_swapped", swaps, "primaries_resolved", resolved);

            // MAPQ rescue and splicing.
            writeMetric(writer, "mapq_zero_in", mMapqZeroIn, "primaries_resolved", resolved);
            writeMetric(writer, "mapq_rescued", mMapqRescued, "mapq_zero_in", mMapqZeroIn);
            writeMetric(writer, "spliced_output", mSplicedOutput, "primaries_resolved", resolved);

            // refinement passes: what each edit actually changed (leading/trailing merged into one figure).
            writeMetric(writer, "supp_merged", mPassEffects.suppMerged());
            writeMetric(writer, "supp_clamped", mPassEffects.suppClamped());
            writeMetric(writer, "overhang_reclaim_records", mPassEffects.reclaimRecords());
            writeMetric(writer, "overhang_reclaim_bases", mPassEffects.reclaimBasesLead() + mPassEffects.reclaimBasesTrail());
            writeMetric(writer, "overhang_collapse", mPassEffects.collapseLeading() + mPassEffects.collapseTrailing());
            writeMetric(writer, "overhang_alt_dropped", mPassEffects.altsDropped());

            // records removed from output.
            writeMetric(writer, "excluded_region_reads", mPassEffects.excludedReads());
            writeMetric(writer, "supp_dropped_low_as", mLowAsSuppsDropped);
            writeMetric(writer, "supp_dropped_orphan", mOrphanSuppsDropped);
            writeMetric(writer, "primary_unmapped_low_as", mLowAsPrimariesUnmapped);
            writeMetric(writer, "primary_unmapped_over_cap", mPassEffects.overCapUnmapped());
        }

        TARS_LOGGER.info("wrote summary to {}", path);
    }

    private static void writeMetric(final BufferedWriter writer, final String metric, final long value) throws IOException
    {
        writer.write(String.join(TSV_DELIM, metric, String.valueOf(value), "", ""));
        writer.newLine();
    }

    private static void writeMetric(
            final BufferedWriter writer, final String metric, final long value, final String basis, final long basisValue)
            throws IOException
    {
        String pct = basisValue > 0 ? String.format("%.2f", 100.0 * value / basisValue) : "";
        writer.write(String.join(TSV_DELIM, metric, String.valueOf(value), pct, basis));
        writer.newLine();
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
                {
                    sb.append(col.name()).append("=").append(count).append(" ");
                }
            }

            if(rowTotal == 0)
                continue;

            TARS_LOGGER.debug("  {}: total={} {}", row.name(), rowTotal, sb.toString());
        }
    }

    static MapqTier deriveMapqTier(final SAMRecord record, final LiftBackResult result)
    {
        int mapq = record.getMappingQuality();
        if(mapq == 0)
        {
            return MapqTier.MAPQ_ZERO;
        }
        return result.numXaAlts() > 0 ? MapqTier.MAPQ_POS_MULTI : MapqTier.MAPQ_POS_UNIQUE;
    }
}
