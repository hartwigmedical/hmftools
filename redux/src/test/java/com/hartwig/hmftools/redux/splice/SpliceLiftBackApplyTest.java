package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
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
                TX_CONTIG, 1, 250, GENE_ID, GENE_NAME, TRANS_NAME, CHR_1, 1,
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

    // Issue #2: willBeUnmapped is the single source of truth for both mate-cache build and the
    // applyAndWriteRecord unmap branches. Cover the cases that previously slipped through.
    @Test
    public void testWillBeUnmappedThresholdLogic()
    {
        LiftBackResolver resolver = new LiftBackResolver(List.of(threeExonContig()));
        SAMRecord r = newRecord(TX_CONTIG, 51, "100M");
        LiftBackResult result = resolver.resolve(r);

        // category UNMAPPED / LIFT_FAILED are always unmapped regardless of threshold
        assertTrue(SpliceLiftBack.willBeUnmapped(
                resolver.resolve(unmappedRecord()), 0, 0));

        // mapped result with thresholds off -> stays mapped
        assertFalse(SpliceLiftBack.willBeUnmapped(result, 0, 0));

        // numLoci=1 < unmap_above_nh=2 -> stays mapped
        assertFalse(SpliceLiftBack.willBeUnmapped(result, 2, 0));

        // numLoci=1 vs unmap_above_nh=0 (off) shouldn't trip
        assertFalse(SpliceLiftBack.willBeUnmapped(result, 0, 0));

        // updatedMapq=60 vs unmap_below_mapq=10 -> stays mapped
        assertFalse(SpliceLiftBack.willBeUnmapped(result, 0, 10));

        // updatedMapq=60 vs unmap_below_mapq=61 -> tripped
        assertTrue(SpliceLiftBack.willBeUnmapped(result, 0, 61));
    }

    @Test
    public void testMarkPrimaryUnmappedClearsAllStaleTags()
    {
        SAMRecord r = newRecord(CHR_1, 1000, "50M");
        r.setReadPairedFlag(true);
        r.setFirstOfPairFlag(true);
        r.setProperPairFlag(true);
        r.setInferredInsertSize(150);
        r.setAttribute("SA", "chrX,1,+,30M,30,0;");
        r.setAttribute("XA", "chrY,+1,50M,0;");
        r.setAttribute("NH", 1);
        r.setAttribute("MC", "50M");
        r.setMappingQuality(60);

        SpliceLiftBack.markPrimaryUnmapped(r);

        assertTrue(r.getReadUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_CIGAR, r.getCigarString());
        assertEquals(0, r.getMappingQuality());
        assertFalse(r.getProperPairFlag());
        assertEquals(0, r.getInferredInsertSize());
        // Every tag that pointed at the now-stale lifted coords is cleared
        assertNull(r.getStringAttribute("SA"));
        assertNull(r.getStringAttribute("XA"));
        assertNull(r.getIntegerAttribute("NH"));
        assertNull(r.getStringAttribute("MC"));
    }

    private static SAMRecord unmappedRecord()
    {
        SAMRecord r = new SAMRecord(new SAMFileHeader());
        r.setReadName("u");
        r.setReadUnmappedFlag(true);
        return r;
    }
}
