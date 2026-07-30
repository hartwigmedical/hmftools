package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

// run-level liftback counters: primaries placed vs lost, and the composition x MAPQ tier cross-tab.
public class LiftBackStatistics
{
    private static final double UNMAPPED_WARN_FRACTION = 0.01;

    public enum MapQualityTier
    {
        MAPQ_ZERO,         // bwa flagged as multi-mapping
        MAPQ_POS_MULTI,    // MAPQ > 0 but XA non-empty (some ambiguity remains)
        MAPQ_POS_UNIQUE    // MAPQ > 0, no XA alts
    }

    private static final int N_COMPOSITIONS = LiftedRecord.Composition.values().length;
    private static final int N_TIERS = MapQualityTier.values().length;

    private final int[][] mCompositionByTier = new int[N_COMPOSITIONS][N_TIERS];
    private int mSwaps = 0;
    private int mTotal = 0;
    private int mResolved = 0;                   // primaries that came out with a genomic placement
    private int mUnmapped = 0;                   // primaries mapped by the aligner that lift-back could not place
    private int mSplicedOutput = 0;              // resolved records emitted with an N (spliced)
    private int mMapQualityZeroIn = 0;           // bwa MAPQ 0 on a resolved primary
    private int mMapQualityRescued = 0;          // of those, lifted to a positive MAPQ

    public void record(final SAMRecord record, final LiftedRecord result)
    {
        // the composition covers the full pre-drop alignment set, so it reflects what bwa produced
        MapQualityTier tier = deriveMapQualityTier(record, result);

        ++mTotal;
        ++mCompositionByTier[result.composition().ordinal()][tier.ordinal()];

        // supplementaries carry a placement too, but the counters below describe primaries
        if(record.getSupplementaryAlignmentFlag())
        {
            return;
        }

        if(!result.hasPlacement())
        {
            // a read that arrived unmapped needed nothing done; only a lost placement is worth counting
            if(!record.getReadUnmappedFlag())
            {
                ++mUnmapped;
            }
            return;
        }

        ++mResolved;

        if(result.swapped())
        {
            ++mSwaps;
        }
        if(result.hasNCigar())
        {
            ++mSplicedOutput;
        }
        if(record.getMappingQuality() == 0)
        {
            ++mMapQualityZeroIn;
            if(result.updatedMapQuality() > 0)
            {
                ++mMapQualityRescued;
            }
        }
    }

    public int total()
    {
        return mTotal;
    }

    public int resolved()
    {
        return mResolved;
    }

    public int unmapped()
    {
        return mUnmapped;
    }

    // folds a worker's counters into this instance for a combined end-of-run summary
    public void merge(final LiftBackStatistics other)
    {
        mTotal += other.mTotal;
        mSplicedOutput += other.mSplicedOutput;
        mMapQualityZeroIn += other.mMapQualityZeroIn;
        mMapQualityRescued += other.mMapQualityRescued;
        mSwaps += other.mSwaps;
        mResolved += other.mResolved;
        mUnmapped += other.mUnmapped;
        for(int i = 0; i < mCompositionByTier.length; ++i)
            for(int j = 0; j < N_TIERS; ++j)
            {
                mCompositionByTier[i][j] += other.mCompositionByTier[i][j];
            }
    }

    public void logSummary()
    {
        TARS_LOGGER.info("processed {} records", mTotal);
        TARS_LOGGER.info("spliced output reads: {}; MAPQ-0 in: {}, rescued: {}", mSplicedOutput, mMapQualityZeroIn, mMapQualityRescued);

        TARS_LOGGER.debug("alignment-set composition x MAPQ tier:");
        logTable(mCompositionByTier, LiftedRecord.Composition.values(), MapQualityTier.values());

        if(mTotal > 0 && mUnmapped / (double) mTotal > UNMAPPED_WARN_FRACTION)
        {
            // ERROR not throw: outputs already written should remain available for inspection
            TARS_LOGGER.error(
                    "unmapped rate {} / {} = {}% exceeds {}% threshold - likely sidecar/FASTA mismatch",
                    mUnmapped, mTotal,
                    String.format("%.2f", 100.0 * mUnmapped / mTotal),
                    String.format("%.2f", 100.0 * UNMAPPED_WARN_FRACTION));
        }
    }

    // flat one-metric-per-line summary: Metric, Value, Pct, Basis. Pct is Value as a percentage of the named
    // Basis metric, so each rate is self-documenting. The cross-tabs stay in the DEBUG log only.
    public void writeSummary(final String path) throws IOException
    {
        try(BufferedWriter writer = createBufferedWriter(path))
        {
            writer.write(String.join(TSV_DELIM, "Metric", "Value", "Pct", "Basis"));
            writer.newLine();

            // primaries_resolved is the denominator for the splits below; a non-zero unmapped rate flags a sidecar/FASTA mismatch
            writeMetric(writer, "records_total", mTotal);
            writeMetric(writer, "primaries_resolved", mResolved, "records_total", mTotal);
            writeMetric(writer, "primaries_unmapped", mUnmapped, "records_total", mTotal);

            writeMetric(writer, "primary_swapped", mSwaps, "primaries_resolved", mResolved);

            writeMetric(writer, "mapQuality_zero_in", mMapQualityZeroIn, "primaries_resolved", mResolved);
            writeMetric(writer, "mapQuality_rescued", mMapQualityRescued, "mapQuality_zero_in", mMapQualityZeroIn);
            writeMetric(writer, "spliced_output", mSplicedOutput, "primaries_resolved", mResolved);
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

    static MapQualityTier deriveMapQualityTier(final SAMRecord record, final LiftedRecord result)
    {
        int mapQuality = record.getMappingQuality();
        if(mapQuality == 0)
        {
            return MapQualityTier.MAPQ_ZERO;
        }
        return result.numXaAlts() > 0 ? MapQualityTier.MAPQ_POS_MULTI : MapQualityTier.MAPQ_POS_UNIQUE;
    }
}
