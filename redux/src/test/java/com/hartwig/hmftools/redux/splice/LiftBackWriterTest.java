package com.hartwig.hmftools.redux.splice;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

// Verifies LiftBackWriter TSV layout: stable headers, one row per record (TSV-A) / per alignment (TSV-B).
public class LiftBackWriterTest
{
    private static SAMRecord pairedRecord(final String readName, final boolean firstOfPair)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(firstOfPair);
        record.setSecondOfPairFlag(!firstOfPair);
        return record;
    }

    private static SAMRecord unpairedRecord(final String readName)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReadPairedFlag(false);
        return record;
    }

    private Path mTsvA;
    private Path mTsvB;

    @Before
    public void setUp() throws IOException
    {
        mTsvA = Files.createTempFile("liftback_a_", ".tsv");
        mTsvB = Files.createTempFile("liftback_b_", ".tsv");
    }

    @After
    public void tearDown() throws IOException
    {
        Files.deleteIfExists(mTsvA);
        Files.deleteIfExists(mTsvB);
    }

    @Test
    public void testHeadersAndRowsForRefOnlyResult() throws IOException
    {
        LiftedAlignment selfAlignment = new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, "1", 1000, "50M",
                "1", 1000, "50M",
                100, 0,
                null, null, null,
                false, true);

        LiftBackResult result = new LiftBackResult(
                LiftBackCategory.REF_SINGLE, LiftBackResult.Composition.REF_ONLY,
                LiftBackResult.RecordRole.PRIMARY,
                "1", 1000, "50M",
                false,
                false, 60, 60,
                0, 1, 0,
                1, 1,
                false, false, false, true,
                "", "",
                0,
                List.of(selfAlignment));

        try(LiftBackWriter writer = new LiftBackWriter(mTsvA.toString(), mTsvB.toString()))
        {
            writer.write(pairedRecord("read1", true), result);
        }

        List<String> aLines = Files.readAllLines(mTsvA);
        List<String> bLines = Files.readAllLines(mTsvB);

        assertEquals(2, aLines.size());
        assertEquals("ReadName\tMateNum\tRecordRole\tPrimaryBucket\tComposition\tDecisionCategory"
                + "\tInputMapq\tUpdatedMapq\tNumXaAlts\tNumRefAlts\tNumTxAlts\tNumLoci"
                + "\tNumDistinctCigarsAtPrimaryLocus\tTxHasNCigar\tTxSoftClipAtBoundary"
                + "\tRefSoftClipped\tRefFullMatch\tFinalChrom\tFinalPos\tFinalCigar\tHasNCigar"
                + "\tGeneIds\tNotes", aLines.get(0));
        assertEquals("read1\t1\tPRIMARY\tRESOLVED\tREF_ONLY\tREF_SINGLE\t60\t60\t0\t1\t0\t1\t1\tfalse\tfalse\tfalse\ttrue\t1\t1000\t50M\tfalse\t\t", aLines.get(1));

        assertEquals(2, bLines.size());
        assertEquals("ReadName\tMateNum\tRecordRole\tSource\tOrigContig\tOrigPos\tOrigCigar"
                + "\tLiftedChrom\tLiftedPos\tLiftedCigar\tAlignmentScore\tNumMismatches"
                + "\tTransName\tGeneId\tGeneName", bLines.get(0));
        assertEquals("read1\t1\tPRIMARY\tSELF\t1\t1000\t50M\t1\t1000\t50M\t100\t0\t\t\t", bLines.get(1));
    }

    @Test
    public void testTsvBHasOneRowPerLiftedAlignment() throws IOException
    {
        LiftedAlignment self = new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, "txContig", 1, "50M",
                "1", 100, "50M",
                100, 0,
                "ENST_X", "ENSG_X", "GENEX",
                false, true);

        LiftedAlignment xa = new LiftedAlignment(
                LiftedAlignment.AlignmentSource.XA_INPUT, "1", 5000, "50M",
                "1", 5000, "50M",
                0, 1,
                null, null, null,
                false, true);

        LiftBackResult result = new LiftBackResult(
                LiftBackCategory.BOTH_MULTI, LiftBackResult.Composition.REF_AND_TX,
                LiftBackResult.RecordRole.PRIMARY,
                "1", 100, "50M",
                false,
                false, 0, 0,
                1, 1, 1,
                2, 1,
                false, false, false, true,
                "ENSG_X", "",
                0,
                List.of(self, xa));

        try(LiftBackWriter writer = new LiftBackWriter(mTsvA.toString(), mTsvB.toString()))
        {
            writer.write(pairedRecord("read2", false), result);
        }

        List<String> bLines = Files.readAllLines(mTsvB);
        assertEquals(3, bLines.size()); // header + 2 alignment rows
        assertTrue(bLines.get(1).startsWith("read2\t2\tPRIMARY\tSELF\t"));
        assertTrue(bLines.get(2).startsWith("read2\t2\tPRIMARY\tXA_INPUT\t"));
    }

    private static LiftBackResult refOnlyResult()
    {
        LiftedAlignment self = new LiftedAlignment(
                LiftedAlignment.AlignmentSource.SELF, "1", 1000, "50M",
                "1", 1000, "50M",
                100, 0,
                null, null, null,
                false, true);

        return new LiftBackResult(
                LiftBackCategory.REF_SINGLE, LiftBackResult.Composition.REF_ONLY,
                LiftBackResult.RecordRole.PRIMARY,
                "1", 1000, "50M",
                false,
                false, 60, 60,
                0, 1, 0,
                1, 1,
                false, false, false, true,
                "", "",
                0,
                List.of(self));
    }

    @Test
    public void headerlessWriterEmitsDataRowsOnly() throws IOException
    {
        try(LiftBackWriter writer = new LiftBackWriter(mTsvA.toString(), mTsvB.toString(), false))
        {
            writer.write(pairedRecord("read1", true), refOnlyResult());
        }

        List<String> aLines = Files.readAllLines(mTsvA);
        List<String> bLines = Files.readAllLines(mTsvB);

        assertEquals(1, aLines.size());
        assertTrue(aLines.get(0).startsWith("read1\t1\tPRIMARY\t"));
        assertEquals(1, bLines.size());
        assertTrue(bLines.get(0).startsWith("read1\t1\tPRIMARY\tSELF\t"));
    }

    @Test
    public void headerLineConstantsMatchWrittenHeaders() throws IOException
    {
        try(LiftBackWriter writer = new LiftBackWriter(mTsvA.toString(), mTsvB.toString()))
        {
            // header-only files
        }

        assertEquals(Files.readAllLines(mTsvA).get(0), LiftBackWriter.TSV_A_HEADER_LINE);
        assertEquals(Files.readAllLines(mTsvB).get(0), LiftBackWriter.TSV_B_HEADER_LINE);
    }

    @Test
    public void concatJoinsShardsUnderSingleHeader() throws IOException
    {
        Path shard0 = Files.createTempFile("shard0_", ".tsv");
        Path shard1 = Files.createTempFile("shard1_", ".tsv");
        Path out = Files.createTempFile("concat_", ".tsv");
        try
        {
            Files.write(shard0, List.of("rowA", "rowB"));
            Files.write(shard1, List.of("rowC"));

            SpliceLiftBack.concatenateTsvShards(
                    List.of(shard0.toString(), shard1.toString()), LiftBackWriter.TSV_A_HEADER_LINE, out.toString());

            List<String> lines = Files.readAllLines(out);
            assertEquals(List.of(LiftBackWriter.TSV_A_HEADER_LINE, "rowA", "rowB", "rowC"), lines);
        }
        finally
        {
            Files.deleteIfExists(shard0);
            Files.deleteIfExists(shard1);
            Files.deleteIfExists(out);
        }
    }

    @Test
    public void testEmptyAlignmentSetStillEmitsTsvARow() throws IOException
    {
        // UNMAPPED: TSV-A row written, TSV-B has header only
        LiftBackResult result = new LiftBackResult(
                LiftBackCategory.UNMAPPED, LiftBackResult.Composition.NONE,
                LiftBackResult.RecordRole.PRIMARY,
                "*", 0, "*",
                false,
                false, 0, 0,
                0, 0, 0,
                0, 0,
                false, false, false, false,
                "", "",
                0,
                List.<LiftedAlignment>of());

        try(LiftBackWriter writer = new LiftBackWriter(mTsvA.toString(), mTsvB.toString()))
        {
            writer.write(unpairedRecord("read3"), result);
        }

        List<String> aLines = Files.readAllLines(mTsvA);
        List<String> bLines = Files.readAllLines(mTsvB);

        assertEquals(2, aLines.size());
        assertTrue(aLines.get(1).startsWith("read3\t0\tPRIMARY\tNON_PRIMARY\tNONE\tUNMAPPED\t"));
        assertEquals(1, bLines.size());
    }
}
