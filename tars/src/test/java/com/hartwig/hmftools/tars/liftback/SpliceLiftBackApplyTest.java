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

// Tests LiftBackRecordOps.applyResultToRecord: chrom/pos/CIGAR rewrite, LIFT_FAILED unmapping, and XA replacement with lifted genomic alts.
public class SpliceLiftBackApplyTest
{
    private static final String XA_TAG = "XA";

    // Local unpaired builder kept on purpose: the shared fixture's primaryRecord sets paired/first-of-pair/proper-pair flags,
    // which these apply-to-record tests must not carry.
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
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnliftableRecordMarkedUnmappedAndStripped()
    {
        // pos 251 is entirely past altEnd (250) - no overhang clamp can save it -> translate fails
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

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
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

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
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertTrue(record.getReadUnmappedFlag());
    }

    // willBeUnmapped is the single source of truth for both mate-cache build and the apply-time unmap branch.
    @Test
    public void testWillBeUnmapped()
    {
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        // an unliftable read is always unmapped
        assertTrue(LiftBackRecordOps.willBeUnmapped(resolver.resolve(unmappedRecord())));

        // clean unique tx read (MAPQ 60, no XA) -> kept
        assertFalse(LiftBackRecordOps.willBeUnmapped(resolver.resolve(newRecord(TX_CONTIG, 51, "100M"))));

        // over the XA cap on a GENOMIC primary: MAPQ 0 with no XA -> unmapped (75+ distinct genomic loci)
        SAMRecord overCap = newRecord(CHR_1, 1000, "100M");
        overCap.setMappingQuality(0);
        assertTrue(LiftBackRecordOps.willBeUnmapped(resolver.resolve(overCap)));

        // a TX-CONTIG primary with MAPQ 0 + no XA hit 75+ transcript contigs of one gene -> one genomic locus,
        // so the REF_ONLY-gated over-cap rule does NOT unmap it; it lifts and is kept.
        SAMRecord txOverCap = newRecord(TX_CONTIG, 51, "100M");
        txOverCap.setMappingQuality(0);
        assertFalse(LiftBackRecordOps.willBeUnmapped(resolver.resolve(txOverCap)));

        // MAPQ 0 but XA present = ordinary few-way multimapper -> kept
        SAMRecord multimapper = newRecord(TX_CONTIG, 51, "100M");
        multimapper.setMappingQuality(0);
        multimapper.setAttribute("XA", CHR_1 + ",+5000,100M,0;");
        assertFalse(LiftBackRecordOps.willBeUnmapped(resolver.resolve(multimapper)));
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

        LiftBackRecordOps.markPrimaryUnmapped(record);

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
        // tx 200 is the last base of exon2: a 1bp exon2 anchor. With ANNOTATED_JUNCTION_MIN_ANCHOR_BP=1
        // the translated junction (1M...N...M) is kept, not rolled into a softclip. This is the resolve+apply
        // layer, before the engine's collapser runs - guards the floor independently of collapse.
        SAMRecord record = newRecord(TX_CONTIG, 200, "50M");
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertLifted(record, CHR_1, 399, "1M100N49M");
    }

    @Test
    public void testTxPrimaryDropsRefXaAltAtSameLocus()
    {
        // tx primary with a softclipped ref XA at the same locus; the discriminator favours tx and drops the ref alt.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testSplicedTxRecordGetsXsAStrand()
    {
        // spliced tx records must carry XS:A:+/- for strand-aware junction interpretation (e.g. Isofox).
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertTrue("expected N in lifted cigar", record.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), record.getAttribute("XS"));
    }

    @Test
    public void testNonSplicedRecordHasNoXsA()
    {
        // a non-spliced genomic record must NOT get XS:A.
        SAMRecord record = newRecord(CHR_1, 1000, "100M");
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertFalse("expected no N in lifted cigar", record.getCigarString().contains("N"));
        assertNull(record.getAttribute("XS"));
    }

    @Test
    public void testStrandSwapReverseComplementsSeqAndQuals()
    {
        // bwa placed the read contiguously on the REVERSE strand (chr1:5000 100M) at MAPQ 0; an XA tx alt crosses a
        // junction on the forward strand and scores higher, so the discriminator swaps to it, flipping the strand
        // reverse->forward. SEQ/quals are stored in the reverse orientation, so applyResultToRecord must
        // reverse-complement the bases and reverse the quals to keep SEQ on the genomic forward strand; otherwise the
        // record is emitted mis-oriented and its recomputed NM inflates to ~read length.
        SAMRecord record = newRecord(CHR_1, 5000, "100M");
        record.setReadNegativeStrandFlag(true);
        record.setMappingQuality(0);
        record.setReadBases(("A".repeat(99) + "C").getBytes());
        byte[] quals = new byte[100];
        for(int i = 0; i < 100; ++i)
            quals[i] = (byte) (i % 40);
        record.setBaseQualities(quals.clone());
        record.setAttribute(XA_TAG, TX_CONTIG + ",+51,100M,0;");

        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        // score the tx alt (index 1) above the ref self (index 0) so the score-based pick swaps to the spliced placement
        LiftBackResult result = resolver.resolve(record, (alns, rec) ->
        {
            alns.get(0).GenomicScore = 10;
            alns.get(1).GenomicScore = 100;
        });
        assertFalse("precondition: swap flips strand to forward", result.negativeStrand());

        LiftBackRecordOps.applyResultToRecord(record, result, new LiftedMateInfoCache());

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
        // no strand change (forward tx primary lifts forward) -> SEQ must be left exactly as-is.
        SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setReadBases(("A".repeat(99) + "C").getBytes());
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

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
