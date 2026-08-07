package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.assertLifted;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BamRecordEmitterTest
{
    private static final String XA_TAG = "XA";

    // The shared fixture sets pair flags, which these emitter tests must not carry.
    private static SAMRecord newRecord(final String contig, final int pos, final String cigar)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReferenceName(contig);
        record.setAlignmentStart(pos);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        return record;
    }

    @Test
    public void testTxPrimaryRewrittenToGenomicCoords()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnliftableRecordMarkedUnmappedAndStripped()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertTrue(record.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, record.getReferenceName());
        assertEquals(SAMRecord.NO_ALIGNMENT_START, record.getAlignmentStart());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, record.getCigarString());
        assertEquals(0, record.getMappingQuality());
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testXaTagRewrittenWithLiftedCoords()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute(XA_TAG, TX_CONTIG + ",+1,50M,0;");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertEquals(CHR_1, record.getReferenceName());
        assertEquals(1000, record.getAlignmentStart());
        assertEquals(CHR_1 + ",+100,50M,0;", record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnmappedRecordLeftUntouched()
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReadUnmappedFlag(true);
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertTrue(record.getReadUnmappedFlag());
    }

    @Test
    public void testWillBeUnmapped()
    {
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        SAMRecord unliftable = unmappedRecord();
        assertTrue(BamRecordEmitter.willBeUnmapped(unliftable, resolver.resolve(unliftable)));

        SAMRecord clean = newRecord(TX_CONTIG, 51, "100M");
        assertFalse(BamRecordEmitter.willBeUnmapped(clean, resolver.resolve(clean)));

        SAMRecord overCap = newRecord(CHR_1, 1000, "100M");
        overCap.setMappingQuality(0);
        assertTrue(BamRecordEmitter.willBeUnmapped(overCap, resolver.resolve(overCap)));

        // Transcript hits may collapse to one genomic locus, so the genomic over-cap rule does not apply.
        SAMRecord txOverCap = newRecord(TX_CONTIG, 51, "100M");
        txOverCap.setMappingQuality(0);
        assertFalse(BamRecordEmitter.willBeUnmapped(txOverCap, resolver.resolve(txOverCap)));

        SAMRecord multimapper = newRecord(TX_CONTIG, 51, "100M");
        multimapper.setMappingQuality(0);
        multimapper.setAttribute("XA", CHR_1 + ",+5000,100M,0;");
        assertFalse(BamRecordEmitter.willBeUnmapped(multimapper, resolver.resolve(multimapper)));
    }

    @Test
    public void testMarkPrimaryUnmappedClearsAllStaleTags()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setProperPairFlag(true);
        record.setInferredInsertSize(150);
        record.setAttribute("SA", "chrX,1,+,30M,30,0;");
        record.setAttribute("XA", "chrY,+1,50M,0;");
        record.setAttribute("NH", 1);
        record.setAttribute("MC", "50M");
        record.setMappingQuality(60);

        BamRecordEmitter.markPrimaryUnmapped(record);

        assertTrue(record.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, record.getCigarString());
        assertEquals(0, record.getMappingQuality());
        assertFalse(record.getProperPairFlag());
        assertEquals(0, record.getInferredInsertSize());
        assertNull(record.getStringAttribute("SA"));
        assertNull(record.getStringAttribute("XA"));
        assertNull(record.getIntegerAttribute("NH"));
        assertNull(record.getStringAttribute("MC"));
    }

    @Test
    public void testTinyAnnotatedJunctionAnchorKept()
    {
        // A genuine 1 bp exon anchor must survive coordinate translation.
        SAMRecord record = newRecord(TX_CONTIG, 200, "50M");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertLifted(record, CHR_1, 399, "1M100N49M");
    }

    @Test
    public void testTxPrimaryDropsRefXaAltAtSameLocus()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testSplicedTxRecordGetsXsAStrand()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertTrue("expected N in lifted cigar", record.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), record.getAttribute("XS"));
    }

    @Test
    public void testNonSplicedRecordHasNoXsA()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "100M");
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertFalse("expected no N in lifted cigar", record.getCigarString().contains("N"));
        assertNull(record.getAttribute("XS"));
    }

    @Test
    public void testStrandSwapReverseComplementsSeqAndQuals()
    {
        // Swapping strand requires SEQ and qualities to follow the emitted genomic orientation.
        SAMRecord record = newRecord(CHR_1, 5000, "100M");
        record.setReadNegativeStrandFlag(true);
        record.setMappingQuality(0);
        record.setReadBases(("A".repeat(99) + "C").getBytes());
        byte[] quals = new byte[100];
        for(int i = 0; i < 100; ++i)
            quals[i] = (byte) (i % 40);
        record.setBaseQualities(quals.clone());
        record.setAttribute(XA_TAG, TX_CONTIG + ",+51,100M,0;");

        LiftedRecord result = TarsTestFixtures.recordBuilder()
                .alignments(List.of(
                        new LiftedAlignment(CHR_1, 5000, "100M", 0, false, false, false, 0),
                        new LiftedAlignment(CHR_1, 150, "50M100N50M", 0, true, false, true, 1)))
                .primaryIndex(1)
                .build();
        assertFalse("precondition: swap flips strand to forward", result.negativeStrand());

        BamRecordEmitter.applyResultToRecord(record, result, new LiftedMatePair(), false);

        assertFalse(record.getReadNegativeStrandFlag());
        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertEquals("G" + "T".repeat(99), record.getReadString());
        byte[] outQuals = record.getBaseQualities();
        assertEquals(quals[99], outQuals[0]);
        assertEquals(quals[0], outQuals[99]);
    }

    @Test
    public void testSameStrandLiftLeavesSeqUntouched()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setReadBases(("A".repeat(99) + "C").getBytes());
        LiftBackDiscriminator resolver = new LiftBackDiscriminator(List.of(threeExonContig()));

        BamRecordEmitter.applyResultToRecord(record, resolver.resolve(record), new LiftedMatePair(), false);

        assertFalse(record.getReadNegativeStrandFlag());
        assertEquals("A".repeat(99) + "C", record.getReadString());
    }

    private static SAMRecord unmappedRecord()
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("u");
        record.setReadUnmappedFlag(true);
        return record;
    }
}
