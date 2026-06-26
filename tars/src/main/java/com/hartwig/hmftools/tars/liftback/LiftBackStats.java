package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.PRIMARY_AS_UNMAP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConstants.SUPP_AS_DROP_THRESHOLD;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

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
            int rescueCandidates, int rescueMerged, int suppClamped,
            int tailEvaluated, int tailExtended, int tailBasesLead, int tailBasesTrail,
            long collapseLeading, long collapseTrailing,
            long junctionsCanonicalized, long overCapUnmapped, long excludedReads)
    {
        public static final PassEffects EMPTY = new PassEffects(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
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
            TARS_LOGGER.info("dropped {} non-rescued supps with source AS < {}",
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

            for(DecidingFeature feature : DecidingFeature.values())
            {
                for(MapqTier tier : MapqTier.values())
                {
                    int count = mFeatureByTier[feature.ordinal()][tier.ordinal()];
                    if(count == 0)
                        continue;
                    writer.write(String.join(TSV_DELIM,
                            "feature_x_mapq", feature.name(), tier.name(), String.valueOf(count)));
                    writer.newLine();
                }
            }

            for(RecordState state : RecordState.values())
            {
                int count = mPerState[state.ordinal()];
                if(count == 0)
                    continue;
                writer.write(String.join(TSV_DELIM, "record_state", state.name(), "COUNT", String.valueOf(count)));
                writer.newLine();
            }

            for(DecidingFeature feature : DecidingFeature.values())
            {
                int swaps = mSwapByFeature[feature.ordinal()];
                if(swaps == 0)
                    continue;
                writer.write(String.join(TSV_DELIM, "feature_x_outcome", feature.name(), "SWAP", String.valueOf(swaps)));
                writer.newLine();
                writer.write(String.join(TSV_DELIM,
                        "feature_x_outcome", feature.name(), "DROP", String.valueOf(mPerFeature[feature.ordinal()] - swaps)));
                writer.newLine();
            }

            for(Outcome outcome : Outcome.values())
            {
                int count = outcomeCount(outcome);
                if(count == 0)
                    continue;
                writeRow(writer, "outcome", outcome.name(), "COUNT", count);
            }

            // headline per-read numbers: how much got placed, spliced, and MAPQ-rescued.
            writeRow(writer, "reads", "total", "COUNT", mTotal);
            writeRow(writer, "reads", "spliced_output", "COUNT", mSplicedOutput);
            writeRow(writer, "reads", "mapq_zero_in", "COUNT", mMapqZeroIn);
            writeRow(writer, "reads", "mapq_rescued", "COUNT", mMapqRescued);

            // per-pass effectiveness - what each refinement pass actually did.
            writeRow(writer, "pass_effect", "rescue", "candidates", mPassEffects.rescueCandidates());
            writeRow(writer, "pass_effect", "rescue", "merged", mPassEffects.rescueMerged());
            writeRow(writer, "pass_effect", "rescue", "supp_clamped", mPassEffects.suppClamped());
            writeRow(writer, "pass_effect", "tail_extend", "evaluated", mPassEffects.tailEvaluated());
            writeRow(writer, "pass_effect", "tail_extend", "extended", mPassEffects.tailExtended());
            writeRow(writer, "pass_effect", "tail_extend", "bases_lead", mPassEffects.tailBasesLead());
            writeRow(writer, "pass_effect", "tail_extend", "bases_trail", mPassEffects.tailBasesTrail());
            writeRow(writer, "pass_effect", "collapse", "leading", mPassEffects.collapseLeading());
            writeRow(writer, "pass_effect", "collapse", "trailing", mPassEffects.collapseTrailing());
            writeRow(writer, "pass_effect", "canonicalize", "shifted", mPassEffects.junctionsCanonicalized());
            writeRow(writer, "pass_effect", "over_cap_unmap", "primaries", mPassEffects.overCapUnmapped());
            writeRow(writer, "pass_effect", "excluded_region", "reads", mPassEffects.excludedReads());
            writeRow(writer, "pass_effect", "low_as", "supps_dropped", mLowAsSuppsDropped);
            writeRow(writer, "pass_effect", "orphan_supp", "dropped", mOrphanSuppsDropped);
            writeRow(writer, "pass_effect", "low_as", "primaries_unmapped", mLowAsPrimariesUnmapped);
        }

        TARS_LOGGER.info("wrote summary to {}", path);
    }

    private static void writeRow(
            final BufferedWriter writer, final String section, final String rowKey, final String columnKey, final long count)
            throws IOException
    {
        writer.write(String.join(TSV_DELIM, section, rowKey, columnKey, String.valueOf(count)));
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
