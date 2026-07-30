package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.tars.liftback.TarsTestFixtures.liftedRecordAt;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// Pair-side lookup semantics (mateOf the opposite side, ownPrimary the same side, null when absent)
// and the mate-field patch those lookups feed: RNEXT/PNEXT/mate-strand/MC and TLEN.
public class LiftedMatePairTest
{
    private static LiftedRecord info(final int alignmentStart)
    {
        return liftedRecordAt(CHR_1, alignmentStart, "50M", false);
    }

    @Test
    public void testMateOfFirstOfPairReturnsSecondInPairInfo()
    {
        LiftedRecord firstInfo = info(100);
        LiftedRecord secondInfo = info(300);
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(true, firstInfo);
        pair.recordPrimary(false, secondInfo);

        assertSame(secondInfo, pair.mateOf(true));
    }

    @Test
    public void testMateOfSecondOfPairReturnsFirstInPairInfo()
    {
        LiftedRecord firstInfo = info(100);
        LiftedRecord secondInfo = info(300);
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(true, firstInfo);
        pair.recordPrimary(false, secondInfo);

        assertSame(firstInfo, pair.mateOf(false));
    }

    @Test
    public void testOwnPrimaryReturnsRecordersOwnSideInfo()
    {
        LiftedRecord firstInfo = info(100);
        LiftedRecord secondInfo = info(300);
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(true, firstInfo);
        pair.recordPrimary(false, secondInfo);

        assertSame(firstInfo, pair.ownPrimary(true));
        assertSame(secondInfo, pair.ownPrimary(false));
    }

    @Test
    public void testNullBeforeAnythingRecorded()
    {
        LiftedMatePair pair = new LiftedMatePair();

        assertNull(pair.mateOf(true));
        assertNull(pair.mateOf(false));
        assertNull(pair.ownPrimary(true));
        assertNull(pair.ownPrimary(false));
    }

    @Test
    public void testNullWhenMateNotYetRecorded()
    {
        // /1 is decided before /2 has been lifted, so its mate lookup is null while its own side is set
        LiftedRecord firstInfo = info(100);
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(true, firstInfo);

        assertNull(pair.mateOf(true));
        assertSame(firstInfo, pair.ownPrimary(true));
    }

    @Test
    public void testPatchSamePairSameChromSetsTlenAndMateFields()
    {
        LiftedMatePair pair = new LiftedMatePair();
        LiftedRecord r2Info = liftedRecordAt("1", 400, "100M", true);
        pair.recordPrimary(false, r2Info);

        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        pair.patchMateFields(r1);

        assertEquals("1", r1.getMateReferenceName());
        assertEquals(400, r1.getMateAlignmentStart());
        assertTrue(r1.getMateNegativeStrandFlag());
        assertFalse(r1.getMateUnmappedFlag());
        // leftmost R1 at 100, R2's 100M ends at 499 -> span = 400, positive on R1
        assertEquals(400, r1.getInferredInsertSize());
        // MC must reflect the mate's lifted CIGAR
        assertEquals("100M", r1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testMateCigarWrittenAsLiftedNCigar()
    {
        // junction-spanning mate: MC must carry the lifted N-CIGAR, not the stale pre-lift value
        LiftedMatePair pair = new LiftedMatePair();
        LiftedRecord r2Info = liftedRecordAt("1", 400, "20M500N30M", true);
        pair.recordPrimary(false, r2Info);

        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        r1.setAttribute(MATE_CIGAR_ATTRIBUTE, "50M"); // stale pre-lift MC
        pair.patchMateFields(r1);

        assertEquals("20M500N30M", r1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testMateCigarClearedWhenMateUnmapped()
    {
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(false, LiftedRecord.unmapped(""));

        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        r1.setAttribute(MATE_CIGAR_ATTRIBUTE, "50M");
        pair.patchMateFields(r1);

        assertNull(r1.getStringAttribute(MATE_CIGAR_ATTRIBUTE));
    }

    @Test
    public void testPatchRightmostReadGetsNegativeTlen()
    {
        LiftedMatePair pair = new LiftedMatePair();
        LiftedRecord r1Info = liftedRecordAt("1", 100, "50M", false);
        pair.recordPrimary(true, r1Info);

        SAMRecord r2 = TarsTestFixtures.pairedRecord("read1", false, "1", 400, "50M", true);
        pair.patchMateFields(r2);

        // span = 449 - 100 + 1 = 350; R2 rightmost -> negative TLEN
        assertEquals(-350, r2.getInferredInsertSize());
    }

    @Test
    public void testPatchDifferentChromsZeroesTlenAndClearsProperPair()
    {
        LiftedMatePair pair = new LiftedMatePair();
        LiftedRecord r2Info = liftedRecordAt("2", 400, "50M", true);
        pair.recordPrimary(false, r2Info);

        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        r1.setProperPairFlag(true);
        pair.patchMateFields(r1);

        assertEquals("2", r1.getMateReferenceName());
        assertEquals(400, r1.getMateAlignmentStart());
        assertEquals(0, r1.getInferredInsertSize());
        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testPatchWithUnmappedMate()
    {
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(false, LiftedRecord.unmapped(""));

        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        r1.setProperPairFlag(true);
        pair.patchMateFields(r1);

        assertTrue(r1.getMateUnmappedFlag());
        // SAM: unmapped mate placed at the read's own position
        assertEquals("1", r1.getMateReferenceName());
        assertEquals(100, r1.getMateAlignmentStart());
        assertEquals(0, r1.getInferredInsertSize());
        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testPatchMissingMateClearsProperPair()
    {
        LiftedMatePair pair = new LiftedMatePair();
        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        r1.setProperPairFlag(true);
        pair.patchMateFields(r1);

        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testUnmappedReadWithMappedMateParkedAtMateLocus()
    {
        // An unmapped read whose mate is mapped is placed at the mate's coordinates so the pair stays
        // together in a coord-sorted BAM. It is not a proper pair and carries no insert size.
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(false, liftedRecordAt("1", 400, "50M", true));

        SAMRecord r1 = TarsTestFixtures.pairedUnmappedRecord("read1", true);
        pair.patchMateFields(r1);

        assertFalse(r1.getMateUnmappedFlag());
        assertEquals("1", r1.getMateReferenceName());
        assertEquals(400, r1.getMateAlignmentStart());
        assertEquals("1", r1.getReferenceName());
        assertEquals(400, r1.getAlignmentStart());
        assertEquals(0, r1.getInferredInsertSize());
        assertFalse(r1.getProperPairFlag());
    }

    @Test
    public void testBothMatesUnmappedClearsOwnCoordinates()
    {
        // When both mates are unmapped, the read's stale pre-lift coordinates are cleared.
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(false, LiftedRecord.unmapped(""));

        SAMRecord r1 = TarsTestFixtures.pairedUnmappedRecord("read1", true);
        r1.setReferenceName("1");          // stale pre-lift placement carried over from bwa
        r1.setAlignmentStart(100);
        pair.patchMateFields(r1);

        assertTrue(r1.getMateUnmappedFlag());
        assertEquals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME, r1.getReferenceName());
        assertEquals(SAMRecord.NO_ALIGNMENT_START, r1.getAlignmentStart());
    }

    @Test
    public void testPatchUnpairedRecordIsNoOp()
    {
        LiftedMatePair pair = new LiftedMatePair();
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read1");
        record.setReadPairedFlag(false);
        record.setReferenceName("1");
        record.setAlignmentStart(100);

        // must not throw or attempt to read first-of-pair flag
        pair.patchMateFields(record);
    }

    @Test
    public void testPatchSuppRecordPatchedFromMatePrimary()
    {
        // supplementary R1 must still use R2's primary info
        LiftedMatePair pair = new LiftedMatePair();
        LiftedRecord r2Info = liftedRecordAt("1", 400, "50M", true);
        pair.recordPrimary(false, r2Info);

        SAMRecord r1Supp = TarsTestFixtures.pairedRecord("read1", true, "1", 200, "30M20S", false);
        r1Supp.setSupplementaryAlignmentFlag(true);
        pair.patchMateFields(r1Supp);

        assertEquals("1", r1Supp.getMateReferenceName());
        assertEquals(400, r1Supp.getMateAlignmentStart());
        assertTrue(r1Supp.getMateNegativeStrandFlag());
    }

    @Test
    public void testPatchOverlappingPairSignFromFivePrime()
    {
        // Overlapping pair at same start: forward 5' (100) < reverse 5' (199), so forward read gets positive TLEN.
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(false, liftedRecordAt("1", 100, "100M", true));

        SAMRecord r1 = TarsTestFixtures.pairedRecord("read1", true, "1", 100, "50M", false);
        pair.patchMateFields(r1);

        assertEquals(100, r1.getInferredInsertSize());
    }

    @Test
    public void testPatchOverlappingPairSecondOfPairForwardIsPositive()
    {
        // Regression: second-of-pair is the forward mate; TLEN sign must follow 5'-end position, not pair order.
        LiftedMatePair pair = new LiftedMatePair();
        pair.recordPrimary(true, liftedRecordAt("1", 100, "39S112M", true));

        SAMRecord r2 = TarsTestFixtures.pairedRecord("read1", false, "1", 100, "117M34S", false);
        pair.patchMateFields(r2);

        // forward 5'=100, reverse 5'=211 -> 211 - 100 + 1 = 112, positive on forward read
        assertEquals(112, r2.getInferredInsertSize());
    }
}
