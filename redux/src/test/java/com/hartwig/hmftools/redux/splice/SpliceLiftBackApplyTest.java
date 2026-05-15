package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import org.junit.Test;

// covers SpliceLiftBack.applyResultToRecord — the SAMRecord-level mutation step. Verifies that
// chrom/pos/CIGAR are rewritten correctly, that LIFT_FAILED marks the record unmapped without
// dropping it, and that XA is replaced with lifted, deduped genomic alts (Stage 4).
public class SpliceLiftBackApplyTest
{
    private static final String GENE_ID = "ENSG_TEST";
    private static final String GENE_NAME = "TESTG";
    private static final String TRANS_NAME = "ENST_TEST";
    private static final String TX_CONTIG = "ens" + GENE_ID + "_" + GENE_NAME + "_" + TRANS_NAME;
    private static final String XA_TAG = "XA";

    private static ContigEntry threeExonContig()
    {
        return new ContigEntry(
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1,
                List.of(new BaseRegion(100, 199), new BaseRegion(300, 399), new BaseRegion(500, 549)));
    }

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

        SpliceLiftBack.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertEquals(CHR_1, record.getReferenceName());
        assertEquals(150, record.getAlignmentStart());
        assertEquals("50M100N50M", record.getCigarString());
        assertNull(record.getStringAttribute(XA_TAG));
    }

    @Test
    public void testUnliftableRecordMarkedUnmappedAndStripped()
    {
        // pos 251 is entirely past altEnd (250) — no overhang clamp can save it -> translate fails
        SAMRecord record = newRecord(TX_CONTIG, 251, "10M");
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));

        SpliceLiftBack.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

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

        SpliceLiftBack.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        // ref primary unchanged; Tx XA alt rewritten to genomic chr1:100/50M
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

        SpliceLiftBack.applyResultToRecord(record, resolver.resolve(record), new LiftedMateInfoCache());

        assertTrue(record.getReadUnmappedFlag());
    }
}
