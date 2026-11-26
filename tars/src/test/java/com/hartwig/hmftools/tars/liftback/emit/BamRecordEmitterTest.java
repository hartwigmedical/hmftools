package com.hartwig.hmftools.tars.liftback.emit;

import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_T0_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULTIMA_TP_TAG;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.EMPTY_REF_GENOME;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.TX_CONTIG;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.assertLifted;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.mappedRecord;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.resolve;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.threeExonContig;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.tars.liftback.EnsemblAnnotationIndex;
import com.hartwig.hmftools.tars.liftback.LiftedAlignment;
import com.hartwig.hmftools.tars.liftback.LiftedRecord;
import com.hartwig.hmftools.tars.liftback.TarsTestFixtures;
import com.hartwig.hmftools.tars.liftback.select.AlignmentSelector;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class BamRecordEmitterTest
{
    private static final String XA_TAG = "XA";

    private static SAMRecord newRecord(final String contig, final int pos, final String cigar)
    {
        return mappedRecord("read", contig, pos, cigar);
    }

    @Test
    public void testTxPrimaryRewrittenToGenomicCoords()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnliftableRecordMarkedUnmappedAndStripped()
    {
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

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
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

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
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertTrue(record.getReadUnmappedFlag());
    }

    @Test
    public void testWillBeUnmapped()
    {
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        SAMRecord unliftable = unmappedRecord();
        assertTrue(BamRecordEmitter.willBeUnmapped(unliftable, resolve(selector, unliftable, EMPTY_REF_GENOME)));

        SAMRecord clean = newRecord(TX_CONTIG, 51, "100M");
        assertFalse(BamRecordEmitter.willBeUnmapped(clean, resolve(selector, clean, EMPTY_REF_GENOME)));

        SAMRecord overCap = newRecord(CHR_1, 1000, "100M");
        overCap.setMappingQuality(0);
        assertTrue(BamRecordEmitter.willBeUnmapped(overCap, resolve(selector, overCap, EMPTY_REF_GENOME)));

        // Transcript hits may collapse to one genomic locus, so the genomic over-cap rule does not apply.
        SAMRecord txOverCap = newRecord(TX_CONTIG, 51, "100M");
        txOverCap.setMappingQuality(0);
        assertFalse(BamRecordEmitter.willBeUnmapped(txOverCap, resolve(selector, txOverCap, EMPTY_REF_GENOME)));

        SAMRecord multimapper = newRecord(TX_CONTIG, 51, "100M");
        multimapper.setMappingQuality(0);
        multimapper.setAttribute("XA", CHR_1 + ",+5000,100M,0;");
        multimapper.setReadBases(TarsTestFixtures.bases("A".repeat(100)));
        RefGenomeInterface ref = new TarsTestFixtures.TestGenome().with(CHR_1, 5200, 'A').asRefGenome();
        assertFalse(BamRecordEmitter.willBeUnmapped(multimapper, resolve(selector, multimapper, ref)));
    }

    @Test
    public void testMarkPrimaryUnmappedUsesReduxSamRepresentation()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);
        record.setReadNegativeStrandFlag(true);
        record.setProperPairFlag(true);
        record.setInferredInsertSize(150);
        record.setAttribute("SA", "chrX,1,+,30M,30,0;");
        record.setAttribute("XA", "chrY,+1,50M,0;");
        record.setAttribute("MC", "50M");
        record.setMappingQuality(60);

        new LiftedMatePair().unmapRead(record);

        assertTrue(record.getReadUnmappedFlag());
        assertTrue("REDUX preserves the alignment strand flag", record.getReadNegativeStrandFlag());
        assertEquals("1:1000", record.getStringAttribute("UM"));
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, record.getCigarString());
        assertEquals(0, record.getMappingQuality());
        assertFalse(record.getProperPairFlag());
        assertEquals(0, record.getInferredInsertSize());
        assertNull(record.getStringAttribute("SA"));
        assertEquals("chrY,+1,50M,0;", record.getStringAttribute("XA"));
        assertNull(record.getStringAttribute("MC"));
    }

    @Test
    public void testTinyAnnotatedJunctionAnchorKept()
    {
        // A genuine 1 bp exon anchor must survive coordinate translation.
        SAMRecord record = newRecord(TX_CONTIG, 200, "50M");
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertLifted(record, CHR_1, 399, "1M100N49M");
    }

    @Test
    public void testTxPrimaryDropsRefXaAltAtSameLocus()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testSplicedTxRecordGetsXsAStrand()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertTrue("expected N in lifted cigar", record.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), record.getAttribute("XS"));
    }

    @Test
    public void testNonSplicedRecordHasNoXsA()
    {
        SAMRecord record = newRecord(CHR_1, 1000, "100M");
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertFalse("expected no N in lifted cigar", record.getCigarString().contains("N"));
        assertNull(record.getAttribute("XS"));
    }

    @Test
    public void testStrandSwapReverseComplementsSeqAndQuals()
    {
        SAMRecord record = newRecord(CHR_1, 5000, "100M");
        record.setReadNegativeStrandFlag(true);
        record.setMappingQuality(0);
        record.setReadBases(("A".repeat(99) + "C").getBytes());
        byte[] quals = new byte[100];
        for(int i = 0; i < 100; ++i)
        {
            quals[i] = (byte) (i % 40);
        }
        record.setBaseQualities(quals.clone());
        record.setAttribute(XA_TAG, TX_CONTIG + ",+51,100M,0;");

        LiftedRecord result = TarsTestFixtures.recordBuilder()
                .alignments(List.of(
                        new LiftedAlignment(CHR_1, 5000, "100M", 0, false, false, 0),
                        new LiftedAlignment(CHR_1, 150, "50M100N50M", 0, true, true, 1)))
                .primaryIndex(1)
                .build();
        assertFalse("precondition: swap flips strand to forward", result.negativeStrand());

        BamRecordEmitter.applyResultToRecord(record, result, new LiftedMatePair());

        assertFalse(record.getReadNegativeStrandFlag());
        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertEquals("G" + "T".repeat(99), record.getReadString());
        byte[] outQuals = record.getBaseQualities();
        assertEquals(quals[99], outQuals[0]);
        assertEquals(quals[0], outQuals[99]);
    }

    @Test
    public void testStrandSwapReversesUltimaTpAndT0()
    {
        SAMRecord record = newRecord(CHR_1, 5000, "100M");
        record.setReadNegativeStrandFlag(true);
        record.setMappingQuality(0);
        byte[] tpValues = new byte[100];
        StringBuilder t0Values = new StringBuilder();
        for(int i = 0; i < 100; ++i)
        {
            tpValues[i] = (byte) (i % 3 - 1);
            t0Values.append((char) ('!' + i % 40));
        }
        record.setAttribute(ULTIMA_TP_TAG, tpValues.clone());
        record.setAttribute(ULTIMA_T0_TAG, t0Values.toString());
        record.setAttribute(XA_TAG, TX_CONTIG + ",+51,100M,0;");

        LiftedRecord result = TarsTestFixtures.recordBuilder()
                .alignments(List.of(
                        new LiftedAlignment(CHR_1, 5000, "100M", 0, false, false, 0),
                        new LiftedAlignment(CHR_1, 150, "50M100N50M", 0, true, true, 1)))
                .primaryIndex(1)
                .build();

        BamRecordEmitter.applyResultToRecord(record, result, new LiftedMatePair());

        assertFalse(record.getReadNegativeStrandFlag());
        assertArrayEquals(reverseArray(tpValues), record.getByteArrayAttribute(ULTIMA_TP_TAG));
        assertEquals(t0Values.reverse().toString(), record.getStringAttribute(ULTIMA_T0_TAG));
    }

    @Test
    public void testSameStrandLiftLeavesSeqUntouched()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setReadBases(("A".repeat(99) + "C").getBytes());
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertFalse(record.getReadNegativeStrandFlag());
        assertEquals("A".repeat(99) + "C", record.getReadString());
    }

    @Test
    public void testWritesRecomputedAlignmentScore()
    {
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("AS", 50);
        LiftedAlignment alignment = new LiftedAlignment(CHR_1, 100, "45M5S", 0, false, true, 0);
        alignment.GenomicScore = 42;
        LiftedRecord result = TarsTestFixtures.recordBuilder().alignments(List.of(alignment)).build();

        BamRecordEmitter.applyResultToRecord(record, result, new LiftedMatePair());

        assertEquals(Integer.valueOf(42), record.getIntegerAttribute("AS"));
    }

    @Test
    public void testClearsStaleAlignmentScoreWhenRescoreIsUnavailable()
    {
        SAMRecord record = newRecord(CHR_1, 100, "50M");
        record.setAttribute("AS", 50);
        LiftedAlignment alignment = new LiftedAlignment(CHR_1, 100, "45M5S", 0, false, true, 0);
        LiftedRecord result = TarsTestFixtures.recordBuilder().alignments(List.of(alignment)).build();

        BamRecordEmitter.applyResultToRecord(record, result, new LiftedMatePair());

        assertNull(record.getAttribute("AS"));
    }

    @Test
    public void testIntronOnlyLiftKeepsInputAlignmentScore()
    {
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute("AS", 100);
        AlignmentSelector selector = new AlignmentSelector(
                List.of(threeExonContig()), EnsemblAnnotationIndex.fromContigEntries(List.of(threeExonContig())));

        BamRecordEmitter.applyResultToRecord(record, resolve(selector, record, EMPTY_REF_GENOME), new LiftedMatePair());

        assertEquals(Integer.valueOf(100), record.getIntegerAttribute("AS"));
    }

    private static SAMRecord unmappedRecord()
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("u");
        record.setReadUnmappedFlag(true);
        return record;
    }
}
