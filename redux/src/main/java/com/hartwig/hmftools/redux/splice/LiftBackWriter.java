package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;

// streaming writer for the two LiftBack TSVs (per-record decision + per-alignment detail).
// owns both BufferedWriters; one instance per SpliceLiftBack run, closed at end.
public class LiftBackWriter implements AutoCloseable
{
    private final BufferedWriter mTsvAWriter;
    private final BufferedWriter mTsvBWriter;

    // MateNum is 1 for first-of-pair, 2 for second-of-pair, 0 for unpaired records.
    // Filtered = "true" when the read pair was dropped from the BAM (e.g. rRNA filter); audit row still emitted.
    private static final String[] TSV_A_HEADERS = {
            "ReadName", "MateNum", "RecordRole", "PrimaryBucket", "Composition", "DecisionCategory",
            "BwaMapq", "UpdatedMapq",
            "NumXaAlts", "NumRefAlts", "NumTxAlts", "NumLoci", "NumDistinctCigarsAtPrimaryLocus",
            "TxHasNCigar", "TxSoftClipAtBoundary", "RefSoftClipped", "RefFullMatch",
            "FinalChrom", "FinalPos", "FinalCigar", "HasNCigar",
            "GeneIds", "Notes", "Filtered"
    };

    private static final String[] TSV_B_HEADERS = {
            "ReadName", "MateNum", "RecordRole", "Source",
            "OrigContig", "OrigPos", "OrigCigar",
            "LiftedChrom", "LiftedPos", "LiftedCigar",
            "AlignmentScore", "NumMismatches",
            "TransName", "GeneId", "GeneName"
    };

    public LiftBackWriter(final String tsvAPath, final String tsvBPath) throws IOException
    {
        mTsvAWriter = createBufferedWriter(tsvAPath);
        mTsvBWriter = createBufferedWriter(tsvBPath);
        writeHeader(mTsvAWriter, TSV_A_HEADERS);
        writeHeader(mTsvBWriter, TSV_B_HEADERS);
    }

    public void write(final SAMRecord record, final LiftBackResult result, final boolean filtered) throws IOException
    {
        final String readName = record.getReadName();
        final String mateNum = mateNumColumn(record);
        writeTsvARow(readName, mateNum, result, filtered);
        for(final LiftedAlignment la : result.LiftedAlignments)
            writeTsvBRow(readName, mateNum, result.Role, la);
    }

    private static String mateNumColumn(final SAMRecord record)
    {
        if(!record.getReadPairedFlag())
            return "0";
        return record.getFirstOfPairFlag() ? "1" : "2";
    }

    private void writeTsvARow(final String readName, final String mateNum, final LiftBackResult result, final boolean filtered)
            throws IOException
    {
        mTsvAWriter.write(String.join(TSV_DELIM,
                readName,
                mateNum,
                result.Role.name(),
                result.Category.primaryBucket().name(),
                result.Comp.name(),
                result.Category.name(),
                String.valueOf(result.BwaMapq),
                String.valueOf(result.UpdatedMapq),
                String.valueOf(result.NumXaAlts),
                String.valueOf(result.NumRefAlts),
                String.valueOf(result.NumTxAlts),
                String.valueOf(result.NumLoci),
                String.valueOf(result.NumDistinctCigarsAtPrimaryLocus),
                String.valueOf(result.TxHasNCigar),
                String.valueOf(result.TxSoftClipAtBoundary),
                String.valueOf(result.RefSoftClipped),
                String.valueOf(result.RefFullMatch),
                result.FinalChrom,
                String.valueOf(result.FinalPos),
                result.FinalCigar,
                String.valueOf(result.HasNCigar),
                nullToEmpty(result.GeneIds),
                nullToEmpty(result.Notes),
                String.valueOf(filtered)));
        mTsvAWriter.newLine();
    }

    private void writeTsvBRow(
            final String readName, final String mateNum, final LiftBackResult.RecordRole role, final LiftedAlignment la)
            throws IOException
    {
        mTsvBWriter.write(String.join(TSV_DELIM,
                readName,
                mateNum,
                role.name(),
                la.Source.name(),
                la.OrigContig,
                String.valueOf(la.OrigPos),
                la.OrigCigar,
                la.LiftedChrom,
                String.valueOf(la.LiftedPos),
                la.LiftedCigar,
                String.valueOf(la.AlignmentScore),
                String.valueOf(la.NumMismatches),
                nullToEmpty(la.TransName),
                nullToEmpty(la.GeneId),
                nullToEmpty(la.GeneName)));
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
