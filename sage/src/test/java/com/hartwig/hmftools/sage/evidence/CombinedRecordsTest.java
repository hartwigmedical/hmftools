package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class CombinedRecordsTest
{
    private static final String REF_BASES = "X" + generateRandomBases(100);

    @Test
    public void testCombinedRecords()
    {
        String readId = "READ_01";
        String chromosome = "1";

        // first a basic match
        SAMRecord first = createSamRecord(readId, chromosome, 1, REF_BASES.substring(1, 21), "20M");

        SAMRecord second = createSamRecord(readId, chromosome, 5, REF_BASES.substring(5, 25), "20M");

        SAMRecord combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(24, combined.getAlignmentEnd());
        assertEquals(REF_BASES.substring(1, 25), combined.getReadString());
        assertEquals("24M", combined.getCigarString());

        // soft-clips extending each end by varying amounts
        first = createSamRecord(readId, chromosome, 6, REF_BASES.substring(1, 21), "5S10M5S");

        second = createSamRecord(readId, chromosome, 12, REF_BASES.substring(10, 30), "2S16M2S");

        combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(6, combined.getAlignmentStart());
        assertEquals(27, combined.getAlignmentEnd());
        assertEquals(REF_BASES.substring(1, 30), combined.getReadString());
        assertEquals("5S22M2S", combined.getCigarString());

        // a delete
        String firstBases = REF_BASES.substring(1, 11) + REF_BASES.substring(12, 22);
        first = createSamRecord(readId, chromosome, 1, firstBases, "10M1D10M");

        String secondBases = REF_BASES.substring(6, 11) + REF_BASES.substring(12, 27);
        second = createSamRecord(readId, chromosome, 6, secondBases, "5M1D15M");

        combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(26, combined.getAlignmentEnd());
        String combinedBases = REF_BASES.substring(1, 11) + REF_BASES.substring(12, 27);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M1D15M", combined.getCigarString());

        // with a longer delete
        firstBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 26);
        first = createSamRecord(readId, chromosome, 1, firstBases, "10M5D10M");

        secondBases = REF_BASES.substring(6, 11) + REF_BASES.substring(16, 31);
        second = createSamRecord(readId, chromosome, 6, secondBases, "5M5D15M");

        combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(30, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 31);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M5D15M", combined.getCigarString());

        // multiple deletes
        firstBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 26) + REF_BASES.substring(28, 38) + REF_BASES.substring(41, 46);
        first = createSamRecord(readId, chromosome, 1, firstBases, "10M5D10M2D10M3D5M");

        secondBases = REF_BASES.substring(6, 11) + REF_BASES.substring(16, 26) + REF_BASES.substring(28, 38) + REF_BASES.substring(41, 51);
        second = createSamRecord(readId, chromosome, 6, secondBases, "5M5D10M2D10M3D10M");

        combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(50, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 26) + REF_BASES.substring(28, 38) + REF_BASES.substring(41, 51);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M5D10M2D10M3D10M", combined.getCigarString());

        // with an insert
        firstBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 21);
        first = createSamRecord(readId, chromosome, 1, firstBases, "10M3I10M");

        secondBases = REF_BASES.substring(6, 11) + "CCC" + REF_BASES.substring(11, 26);
        second = createSamRecord(readId, chromosome, 6, secondBases, "5M3I15M");

        combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(25, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 26);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M3I15M", combined.getCigarString());

        // more complicated example
        firstBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 21) + REF_BASES.substring(26, 36) + "AA"
                + REF_BASES.substring(36, 46);

        first = createSamRecord(readId, chromosome, 1, firstBases, "10M3I10M5D10M2I5M5S");

        secondBases = REF_BASES.substring(6, 11) + "CCC" + REF_BASES.substring(11, 21) + REF_BASES.substring(26, 36) + "AA"
                + REF_BASES.substring(36, 51);
        second = createSamRecord(readId, chromosome, 6, secondBases, "5M3I10M5D10M2I15M");

        combined = formFragmentRead(first, second);
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(50, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 21) + REF_BASES.substring(26, 36) + "AA"
                + REF_BASES.substring(36, 51);

        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M3I10M5D10M2I15M", combined.getCigarString());
    }

    @Test
    public void testMismatches()
    {
        String readId = "READ_01";
        String chromosome = "1";

        // first a basic match
        SAMRecord first = createSamRecord(
                readId, chromosome, 1, REF_BASES.substring(1, 11) + "C" + REF_BASES.substring(11, 21), "10M1I10M");

        SAMRecord second = createSamRecord(readId, chromosome, 5, REF_BASES.substring(5, 25), "20M");

        SyncFragmentOutcome syncOutcome = ReadContextEvidence.formFragmentRead(first, second);
        assertEquals(SyncFragmentType.CIGAR_MISMATCH, syncOutcome.SyncType);

        // off by 1
        first = createSamRecord(
                readId, chromosome, 1, REF_BASES.substring(1, 11) + "C" + REF_BASES.substring(11, 21), "10M1I10M");

        second = createSamRecord(
                readId, chromosome, 2, REF_BASES.substring(1, 12) + "C" + REF_BASES.substring(12, 21), "10M1I10M");

        syncOutcome = ReadContextEvidence.formFragmentRead(first, second);
        assertEquals(SyncFragmentType.CIGAR_MISMATCH, syncOutcome.SyncType);

        // too many mismatches
        first = createSamRecord(readId, chromosome, 1, REF_BASES.substring(1, 21), "20M");
        second = createSamRecord(readId, chromosome, 1, REF_BASES.substring(2, 22), "20M");

        syncOutcome = ReadContextEvidence.formFragmentRead(first, second);
        assertEquals(SyncFragmentType.BASE_MISMATCH, syncOutcome.SyncType);

        // non-overlapping but different INDELs
        first = createSamRecord(
                readId, chromosome, 1, REF_BASES.substring(1, 6) + "C" + REF_BASES.substring(6, 41), "5M1I35M");

        second = createSamRecord(
                readId, chromosome, 30, REF_BASES.substring(30, 40) + REF_BASES.substring(45, 75), "10M5D30M");

        syncOutcome = ReadContextEvidence.formFragmentRead(first, second);
        assertEquals(SyncFragmentType.NO_OVERLAP_CIGAR_DIFF, syncOutcome.SyncType);
    }

        @Test
    public void testCompatibleCigars()
    {
        Cigar first = new Cigar();
        first.add(new CigarElement(10, S));
        first.add(new CigarElement(30, M));
        first.add(new CigarElement(2, I));
        first.add(new CigarElement(20, M));
        first.add(new CigarElement(5, D));
        first.add(new CigarElement(40, M));
        first.add(new CigarElement(8, S));

        assertTrue(ReadContextEvidence.compatibleCigars(first, first));

        // other diffs are not permitted
        Cigar second = new Cigar();
        second.add(new CigarElement(30, M));
        second.add(new CigarElement(3, I));
        second.add(new CigarElement(20, M));
        second.add(new CigarElement(5, D));
        second.add(new CigarElement(40, M));

        assertFalse(ReadContextEvidence.compatibleCigars(first, second));

        second.add(new CigarElement(30, M));
        second.add(new CigarElement(13, D));
        second.add(new CigarElement(40, M));

        assertFalse(ReadContextEvidence.compatibleCigars(first, second));

        // can differ in soft-clips and aligned lengths
        first = new Cigar();
        first.add(new CigarElement(10, S));
        first.add(new CigarElement(30, M));
        first.add(new CigarElement(12, S));

        second = new Cigar();
        second.add(new CigarElement(8, S));
        second.add(new CigarElement(40, M));
        second.add(new CigarElement(2, S));

        assertTrue(ReadContextEvidence.compatibleCigars(first, first));

        second = new Cigar();
        second.add(new CigarElement(40, M));

        assertTrue(ReadContextEvidence.compatibleCigars(first, first));
    }

    private static SAMRecord formFragmentRead(final SAMRecord first, final SAMRecord second)
    {
        return ReadContextEvidence.formFragmentRead(first, second).CombinedRecord;
    }

}
