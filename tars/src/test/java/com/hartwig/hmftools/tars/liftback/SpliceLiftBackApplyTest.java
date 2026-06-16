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
        // pos 251 is entirely past altEnd (250) — no overhang clamp can save it -> translate fails
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

    // Issue #2: willBeUnmapped is the single source of truth for both mate-cache build and the
    // applyAndWriteRecord unmap branches. Cover the cases that previously slipped through.
    @Test
    public void testWillBeUnmappedThresholdLogic()
    {
        final LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        final SAMRecord record = newRecord(TX_CONTIG, 51, "100M");
        final LiftBackResult result = resolver.resolve(record);

        assertTrue(LiftBackRecordOps.willBeUnmapped(
                resolver.resolve(unmappedRecord()), 0, 0));

        assertFalse(LiftBackRecordOps.willBeUnmapped(result, 0, 0));

        assertFalse(LiftBackRecordOps.willBeUnmapped(result, 2, 0));

        assertFalse(LiftBackRecordOps.willBeUnmapped(result, 0, 0));

        assertFalse(LiftBackRecordOps.willBeUnmapped(result, 0, 10));

        assertTrue(LiftBackRecordOps.willBeUnmapped(result, 0, 61));
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

    private static SAMRecord unmappedRecord()
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("u");
        record.setReadUnmappedFlag(true);
        return record;
    }
}
