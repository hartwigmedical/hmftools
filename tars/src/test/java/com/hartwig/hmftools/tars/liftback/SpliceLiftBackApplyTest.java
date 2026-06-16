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
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
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
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnliftableRecordMarkedUnmappedAndStripped()
    {
        // pos 251 is entirely past altEnd (250) - no overhang clamp can save it -> translate fails
        final SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

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
        final SAMRecord record = newRecord(CHR_1, 1000, "50M");
        record.setAttribute(XA_TAG, TX_CONTIG + ",+1,50M,0;");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertEquals(CHR_1, record.getReferenceName());
        assertEquals(1000, record.getAlignmentStart());
        assertEquals(CHR_1 + ",+100,50M,0;", record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnmappedRecordLeftUntouched()
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReadUnmappedFlag(true);
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertTrue(record.getReadUnmappedFlag());
    }

    // willBeUnmapped is the single source of truth for both mate-cache build and the apply-time unmap branch.
    @Test
    public void testWillBeUnmapped()
    {
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        // an unliftable read is always unmapped
        assertTrue(LiftBackRecordOps.willBeUnmapped(resolver.resolve(unmappedRecord())));

        // clean unique tx read (MAPQ 60, no XA) -> kept
        assertFalse(LiftBackRecordOps.willBeUnmapped(resolver.resolve(newRecord(TX_CONTIG, 51, "100M"))));

        // over the XA cap: MAPQ 0 with no XA -> unmapped (maps to too many loci to place)
        final SAMRecord overCap = newRecord(TX_CONTIG, 51, "100M");
        overCap.setMappingQuality(0);
        assertTrue(LiftBackRecordOps.willBeUnmapped(resolver.resolve(overCap)));

        // MAPQ 0 but XA present = ordinary few-way multimapper -> kept
        final SAMRecord multimapper = newRecord(TX_CONTIG, 51, "100M");
        multimapper.setMappingQuality(0);
        multimapper.setAttribute("XA", CHR_1 + ",+5000,100M,0;");
        assertFalse(LiftBackRecordOps.willBeUnmapped(resolver.resolve(multimapper)));
    }

    @Test
    public void testMarkPrimaryUnmappedClearsAllStaleTags()
    {
        final SAMRecord record = newRecord(CHR_1, 1000, "50M");
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
        final SAMRecord record = newRecord(TX_CONTIG, 200, "50M");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertLifted(record, CHR_1, 399, "1M100N49M");
    }

    @Test
    public void testTxPrimaryDropsRefXaAltAtSameLocus()
    {
        // tx primary with a softclipped ref XA at the same locus; the discriminator favours tx and drops the ref alt.
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        record.setAttribute(XA_TAG, CHR_1 + ",+150,50M50S,0;");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertLifted(record, CHR_1, 150, "50M100N50M");
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testSplicedTxRecordGetsXsAStrand()
    {
        // spliced tx records must carry XS:A:+/- for strand-aware junction interpretation (e.g. Isofox).
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertTrue("expected N in lifted cigar", record.getCigarString().contains("N"));
        assertEquals(Character.valueOf('+'), record.getAttribute("XS"));
    }

    @Test
    public void testNonSplicedRecordHasNoXsA()
    {
        // a non-spliced genomic record must NOT get XS:A.
        final SAMRecord record = newRecord(CHR_1, 1000, "100M");
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        LiftBackRecordOps.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertFalse("expected no N in lifted cigar", record.getCigarString().contains("N"));
        assertNull(record.getAttribute("XS"));
    }

    private static SAMRecord unmappedRecord()
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("u");
        record.setReadUnmappedFlag(true);
        return record;
    }
}
