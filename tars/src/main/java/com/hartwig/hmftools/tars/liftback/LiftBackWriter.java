package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

public class LiftBackWriter implements AutoCloseable
{
    private final BufferedWriter mTsvAWriter;
    private final BufferedWriter mTsvBWriter;

    private static final String[] TSV_A_HEADERS = {
            "ReadName", "MateNum", "RecordRole", "RecordState", "Composition", "Outcome", "DecidingFeature", "Swapped",
            "InputMapq", "UpdatedMapq",
            "NumXaAlts", "NumRefAlts", "NumTxAlts", "NumLoci", "NumDistinctCigarsAtPrimaryLocus",
            "TxHasNCigar", "TxSoftClipAtBoundary", "RefSoftClipped", "RefFullMatch",
            "FinalChrom", "FinalPos", "FinalCigar", "HasNCigar",
            "GeneIds", "Notes"
    };

    private static final String[] TSV_B_HEADERS = {
            "ReadName", "MateNum", "RecordRole", "Source",
            "OrigContig", "OrigPos", "OrigCigar",
            "LiftedChrom", "LiftedPos", "LiftedCigar",
            "AlignmentScore", "NumMismatches",
            "TransName", "GeneId", "GeneName"
    };

    // header lines for the records (A) and alignments (B) TSVs, written once by the shard-concat step.
    public static final String TSV_A_HEADER_LINE = String.join(TSV_DELIM, TSV_A_HEADERS);
    public static final String TSV_B_HEADER_LINE = String.join(TSV_DELIM, TSV_B_HEADERS);

    public LiftBackWriter(final String tsvAPath, final String tsvBPath) throws IOException
    {
        this(tsvAPath, tsvBPath, true);
    }

    // writeHeaders=false produces a data-only shard: each worker writes one, and the driver concatenates the
    // shards under a single header line (see TarsApplication.concatenateTsvShards).
    public LiftBackWriter(final String tsvAPath, final String tsvBPath, final boolean writeHeaders) throws IOException
    {
        mTsvAWriter = createBufferedWriter(tsvAPath);
        mTsvBWriter = createBufferedWriter(tsvBPath);
        if(writeHeaders)
        {
            writeHeader(mTsvAWriter, TSV_A_HEADERS);
            writeHeader(mTsvBWriter, TSV_B_HEADERS);
        }
    }

    public void write(final SAMRecord record, final LiftBackResult result) throws IOException
    {
        String readName = record.getReadName();
        String mateNum = mateNumColumn(record);
        writeTsvARow(readName, mateNum, result);
        for(LiftedAlignment alignment : result.liftedAlignments())
        {
            writeTsvBRow(readName, mateNum, result.role(), alignment);
        }
    }

    private static String mateNumColumn(final SAMRecord record)
    {
        if(!record.getReadPairedFlag())
        {
            return "0";
        }
        return record.getFirstOfPairFlag() ? "1" : "2";
    }

    private void writeTsvARow(final String readName, final String mateNum, final LiftBackResult result) throws IOException
    {
        mTsvAWriter.write(String.join(TSV_DELIM,
                readName,
                mateNum,
                result.role().name(),
                result.recordState().name(),
                result.comp().name(),
                result.outcome().name(),
                result.decidingFeature() != null ? result.decidingFeature().name() : "",
                String.valueOf(result.swapped()),
                String.valueOf(result.inputMapq()),
                String.valueOf(result.updatedMapq()),
                String.valueOf(result.numXaAlts()),
                String.valueOf(result.numRefAlts()),
                String.valueOf(result.numTxAlts()),
                String.valueOf(result.numLoci()),
                String.valueOf(result.numDistinctCigarsAtPrimaryLocus()),
                String.valueOf(result.txHasNCigar()),
                String.valueOf(result.txSoftClipAtBoundary()),
                String.valueOf(result.refSoftClipped()),
                String.valueOf(result.refFullMatch()),
                result.finalChrom(),
                String.valueOf(result.finalPos()),
                result.finalCigar(),
                String.valueOf(result.hasNCigar()),
                nullToEmpty(result.geneIds()),
                nullToEmpty(result.notes())));
        mTsvAWriter.newLine();
    }

    private void writeTsvBRow(
            final String readName, final String mateNum, final LiftBackResult.RecordRole role, final LiftedAlignment alignment)
            throws IOException
    {
        mTsvBWriter.write(String.join(TSV_DELIM,
                readName,
                mateNum,
                role.name(),
                alignment.Source.name(),
                alignment.OrigContig,
                String.valueOf(alignment.OrigPos),
                alignment.OrigCigar,
                alignment.LiftedChrom,
                String.valueOf(alignment.LiftedPos),
                alignment.LiftedCigar,
                String.valueOf(alignment.AlignmentScore),
                String.valueOf(alignment.NumMismatches),
                nullToEmpty(alignment.TransName),
                nullToEmpty(alignment.GeneId),
                nullToEmpty(alignment.GeneName)));
        mTsvBWriter.newLine();
    }

    private static void writeHeader(final BufferedWriter writer, final String[] headers) throws IOException
    {
        writer.write(String.join(TSV_DELIM, headers));
        writer.newLine();
    }

    private static String nullToEmpty(final String value)
    {
        return value != null ? value : "";
    }

    @Override
    public void close() throws IOException
    {
        mTsvAWriter.close();
        mTsvBWriter.close();
    }
}
