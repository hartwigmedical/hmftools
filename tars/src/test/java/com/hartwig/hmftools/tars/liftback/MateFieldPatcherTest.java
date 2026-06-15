package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class MateFieldPatcherTest
{
    private static SAMRecord pairedMappedRecord(
            final String readName, final boolean firstOfPair,
            final String chromosome, final int alignmentStart, final String cigar, final boolean negStrand)
    {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(readName);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(firstOfPair);
        record.setSecondOfPairFlag(!firstOfPair);
        record.setReadUnmappedFlag(false);
        record.setReferenceName(chromosome);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadNegativeStrandFlag(negStrand);
        record.setMappingQuality(60);
        return record;
    }

    @Test
    public void testPatchSamePairSameChromSetsTlenAndMateFields()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        LiftedMateInfo r2Info = LiftedMateInfo.mapped("1", 400, 499, "50M", true);
        cache.recordPrimaryAlignment("read1", false, r2Info);

        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        MateFieldPatcher.patchMateFields(r1, cache);

        assertEquals("1", r1.getMateReferenceName());
        assertEquals(400, r1.getMateAlignmentStart());
        assertTrue(r1.getMateNegativeStrandFlag());
        assertFalse(r1.getMateUnmappedFlag());
        // leftmost R1 at 100, rightmost end R2 at 499 → span = 400, positive on R1
        assertEquals(400, r1.getInferredInsertSize());
        // MC must reflect the partner's lifted CIGAR
        assertEquals("50M", r1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testMateCigarWrittenAsLiftedNCigar()
    {
        // junction-spanning partner: MC must carry the lifted N-CIGAR, not the stale pre-lift value
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        LiftedMateInfo r2Info = LiftedMateInfo.mapped("1", 400, 999, "20M500N30M", true);
        cache.recordPrimaryAlignment("read1", false, r2Info);

        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        r1.setAttribute(MATE_CIGAR_ATTRIBUTE, "50M"); // stale pre-lift MC
        MateFieldPatcher.patchMateFields(r1, cache);

        assertEquals("20M500N30M", r1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testMateCigarClearedWhenPartnerUnmapped()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        cache.recordPrimaryAlignment("read1", false, LiftedMateInfo.UNMAPPED);

        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        r1.setAttribute(MATE_CIGAR_ATTRIBUTE, "50M");
        MateFieldPatcher.patchMateFields(r1, cache);

        assertNull(r1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testPatchRightmostReadGetsNegativeTlen()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        LiftedMateInfo r1Info = LiftedMateInfo.mapped("1", 100, 149, "50M", false);
        cache.recordPrimaryAlignment("read1", true, r1Info);

        SAMRecord r2 = pairedMappedRecord("read1", false, "1", 400, "50M", true);
        MateFieldPatcher.patchMateFields(r2, cache);

        // span = 449 - 100 + 1 = 350; R2 rightmost → negative TLEN
        assertEquals(-350, r2.getInferredInsertSize());
    }

    @Test
    public void testPatchDifferentChromsZeroesTlenAndClearsProperPair()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        LiftedMateInfo r2Info = LiftedMateInfo.mapped("2", 400, 499, "50M", true);
        cache.recordPrimaryAlignment("read1", false, r2Info);

        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        r1.setProperPairFlag(true);
        MateFieldPatcher.patchMateFields(r1, cache);

        assertEquals("2", r1.getMateReferenceName());
        assertEquals(400, r1.getMateAlignmentStart());
        assertEquals(0, r1.getInferredInsertSize());
        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testPatchWithUnmappedPartner()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        cache.recordPrimaryAlignment("read1", false, LiftedMateInfo.UNMAPPED);

        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        r1.setProperPairFlag(true);
        MateFieldPatcher.patchMateFields(r1, cache);

        assertTrue(r1.getMateUnmappedFlag());
        // SAM: unmapped mate placed at the read's own position
        assertEquals("1", r1.getMateReferenceName());
        assertEquals(100, r1.getMateAlignmentStart());
        assertEquals(0, r1.getInferredInsertSize());
        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testPatchMissingPartnerClearsProperPair()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        r1.setProperPairFlag(true);
        MateFieldPatcher.patchMateFields(r1, cache);

        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testPatchUnpairedRecordIsNoOp()
    {
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read1");
        record.setReadPairedFlag(false);
        record.setReferenceName("1");
        record.setAlignmentStart(100);

        // must not throw or attempt to read first-of-pair flag
        MateFieldPatcher.patchMateFields(record, cache);
    }

    @Test
    public void testPatchSuppRecordPatchedFromPartnerPrimary()
    {
        // supplementary R1 must still use R2's primary info
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        LiftedMateInfo r2Info = LiftedMateInfo.mapped("1", 400, 499, "50M", true);
        cache.recordPrimaryAlignment("read1", false, r2Info);

        SAMRecord r1Supp = pairedMappedRecord("read1", true, "1", 200, "30M20S", false);
        r1Supp.setSupplementaryAlignmentFlag(true);
        MateFieldPatcher.patchMateFields(r1Supp, cache);

        assertEquals("1", r1Supp.getMateReferenceName());
        assertEquals(400, r1Supp.getMateAlignmentStart());
        assertTrue(r1Supp.getMateNegativeStrandFlag());
    }

    @Test
    public void testPatchOverlappingPairSignFromFivePrime()
    {
        // Overlapping pair at same start: forward 5' (100) < reverse 5' (199), so forward read gets positive TLEN.
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        cache.recordPrimaryAlignment("read1", false, LiftedMateInfo.mapped("1", 100, 199, "50M", true));

        SAMRecord r1 = pairedMappedRecord("read1", true, "1", 100, "50M", false);
        MateFieldPatcher.patchMateFields(r1, cache);

        assertEquals(100, r1.getInferredInsertSize());
    }

    @Test
    public void testPatchOverlappingPairSecondOfPairForwardIsPositive()
    {
        // Regression: second-of-pair is the forward mate; TLEN sign must follow 5'-end position, not pair order.
        LiftedMateInfoCache cache = new LiftedMateInfoCache();
        cache.recordPrimaryAlignment("read1", true, LiftedMateInfo.mapped("1", 100, 211, "39S112M", true));

        SAMRecord r2 = pairedMappedRecord("read1", false, "1", 100, "117M34S", false);
        MateFieldPatcher.patchMateFields(r2, cache);

        // forward 5'=100, reverse 5'=211 → 211 - 100 + 1 = 112, positive on forward read
        assertEquals(112, r2.getInferredInsertSize());
    }
}
