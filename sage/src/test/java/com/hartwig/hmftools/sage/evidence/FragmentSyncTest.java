package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.test.MockRefGenome.getNextBase;
import static com.hartwig.hmftools.sage.common.TestUtils.MOCK_REF_GENOME;
import static com.hartwig.hmftools.sage.common.TestUtils.MSI_JITTER_CALCS;
import static com.hartwig.hmftools.sage.common.TestUtils.RECALIBRATION;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.sync.CombinedSyncData.formFragmentRead;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.sync.CombinedSyncData;
import com.hartwig.hmftools.sage.sync.FragmentData;
import com.hartwig.hmftools.sage.sync.FragmentSyncOutcome;
import com.hartwig.hmftools.sage.sync.FragmentSyncType;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class FragmentSyncTest
{
    private static final String REF_BASES = "X" + generateRandomBases(100);
    private static final String READ_ID = "READ_01";

    @Test
    public void testCombinedRecords()
    {
        // first a basic match
        SAMRecord first = createSamRecord(READ_ID, CHR_1, 1, REF_BASES.substring(1, 21), "20M");

        SAMRecord second = createSamRecord(READ_ID, CHR_1, 5, REF_BASES.substring(5, 25), "20M");
        second.setReadNegativeStrandFlag(true);

        SAMRecord combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(24, combined.getAlignmentEnd());
        assertEquals(REF_BASES.substring(1, 25), combined.getReadString());
        assertEquals("24M", combined.getCigarString());

        // soft-clips extending each end by varying amounts
        first = createSamRecord(READ_ID, CHR_1, 6, REF_BASES.substring(1, 21), "5S10M5S");

        second = createSamRecord(READ_ID, CHR_1, 12, REF_BASES.substring(10, 30), "2S16M2S");
        second.setReadNegativeStrandFlag(true);

        combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(6, combined.getAlignmentStart());
        assertEquals(27, combined.getAlignmentEnd());
        assertEquals(REF_BASES.substring(1, 30), combined.getReadString());
        assertEquals("5S22M2S", combined.getCigarString());

        // a delete
        String firstBases = REF_BASES.substring(1, 11) + REF_BASES.substring(12, 22);
        first = createSamRecord(READ_ID, CHR_1, 1, firstBases, "10M1D10M");

        String secondBases = REF_BASES.substring(6, 11) + REF_BASES.substring(12, 27);
        second = createSamRecord(READ_ID, CHR_1, 6, secondBases, "5M1D15M");
        second.setReadNegativeStrandFlag(true);

        combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(26, combined.getAlignmentEnd());
        String combinedBases = REF_BASES.substring(1, 11) + REF_BASES.substring(12, 27);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M1D15M", combined.getCigarString());

        // with a longer delete
        firstBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 26);
        first = createSamRecord(READ_ID, CHR_1, 1, firstBases, "10M5D10M");

        secondBases = REF_BASES.substring(6, 11) + REF_BASES.substring(16, 31);
        second = createSamRecord(READ_ID, CHR_1, 6, secondBases, "5M5D15M");
        second.setReadNegativeStrandFlag(true);

        combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(30, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 31);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M5D15M", combined.getCigarString());

        // multiple deletes
        firstBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 26) + REF_BASES.substring(28, 38) + REF_BASES.substring(41, 46);
        first = createSamRecord(READ_ID, CHR_1, 1, firstBases, "10M5D10M2D10M3D5M");

        secondBases = REF_BASES.substring(6, 11) + REF_BASES.substring(16, 26) + REF_BASES.substring(28, 38) + REF_BASES.substring(41, 51);
        second = createSamRecord(READ_ID, CHR_1, 6, secondBases, "5M5D10M2D10M3D10M");
        second.setReadNegativeStrandFlag(true);

        combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(50, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + REF_BASES.substring(16, 26) + REF_BASES.substring(28, 38) + REF_BASES.substring(41, 51);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M5D10M2D10M3D10M", combined.getCigarString());

        // with an insert
        firstBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 21);
        first = createSamRecord(READ_ID, CHR_1, 1, firstBases, "10M3I10M");

        secondBases = REF_BASES.substring(6, 11) + "CCC" + REF_BASES.substring(11, 26);
        second = createSamRecord(READ_ID, CHR_1, 6, secondBases, "5M3I15M");
        second.setReadNegativeStrandFlag(true);

        combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(25, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 26);
        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M3I15M", combined.getCigarString());

        // more complicated example
        firstBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 21) + REF_BASES.substring(26, 36) + "AA"
                + REF_BASES.substring(36, 46);

        first = createSamRecord(READ_ID, CHR_1, 1, firstBases, "10M3I10M5D10M2I5M5S");

        secondBases = REF_BASES.substring(6, 11) + "CCC" + REF_BASES.substring(11, 21) + REF_BASES.substring(26, 36) + "AA"
                + REF_BASES.substring(36, 51);
        second = createSamRecord(READ_ID, CHR_1, 6, secondBases, "5M3I10M5D10M2I15M");
        second.setReadNegativeStrandFlag(true);

        combined = formFragmentRead(first, second).CombinedRecord;
        assertNotNull(combined);
        assertEquals(1, combined.getAlignmentStart());
        assertEquals(50, combined.getAlignmentEnd());
        combinedBases = REF_BASES.substring(1, 11) + "CCC" + REF_BASES.substring(11, 21) + REF_BASES.substring(26, 36) + "AA"
                + REF_BASES.substring(36, 51);

        assertEquals(combinedBases, combined.getReadString());
        assertEquals("10M3I10M5D10M2I15M", combined.getCigarString());
    }

    @Test
    public void testCombinedRecords2()
    {
        SAMRecord first = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 25) + REF_BASES.substring(35, 55), "15M10D20M");

        SAMRecord second = createSamRecord(READ_ID, CHR_1, 45, REF_BASES.substring(45, 75), "20M10S");
        second.setReadNegativeStrandFlag(true);

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        SAMRecord combined = syncOutcome.CombinedRecord;
        assertNotNull(combined);
        assertEquals(10, combined.getAlignmentStart());
        assertEquals(64, combined.getAlignmentEnd());
        assertEquals("15M10D30M10S", combined.getCigarString());
        assertEquals(REF_BASES.substring(10, 25) + REF_BASES.substring(35, 75), combined.getReadString());
    }

    @Test
    public void testSoftClips()
    {
        SAMRecord first = createSamRecord(
                READ_ID, CHR_1, 20, REF_BASES.substring(20, 35), "10M5S");
        first.setReadNegativeStrandFlag(true);

        SAMRecord second = createSamRecord(
                READ_ID, CHR_1, 20, REF_BASES.substring(17, 30), "3S10M");

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        SAMRecord combined = syncOutcome.CombinedRecord;
        assertNotNull(combined);
        assertEquals(20, combined.getAlignmentStart());
        assertEquals(29, combined.getAlignmentEnd());
        assertEquals("3S10M5S", combined.getCigarString());
        assertEquals(REF_BASES.substring(17, 35), combined.getReadString());
    }

    @Test
    public void testReadStrandMismatches()
    {
        // one read supports the variant, one doesn't and this factors into the read strand determination
        int position = 20;

        String refBase = REF_BASES.substring(position, position + 1);
        String altBase = getNextBase(refBase);

        SimpleVariant variant = createSimpleVariant(position, refBase, altBase);

        VariantReadContext readContext = createReadContext(
                variant, REF_BASES.substring(position - 2, position), REF_BASES.substring(position + 1, position + 2));

        final RefSequence refSequence = new RefSequence(1, REF_BASES.getBytes());

        QualityCalculator qualityCalculator = new QualityCalculator(
                TEST_CONFIG, RECALIBRATION, refSequence, MOCK_REF_GENOME, MSI_JITTER_CALCS);

        ReadContextCounter readContextCounter = createReadCounter(readContext);

        String READ_ID = "READ_01";

        SAMRecord first = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 40), "30M");
        first.getBaseQualities()[10] = (byte)11;

        SAMRecord second = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, position) + altBase + REF_BASES.substring(position + 1, 40), "30M");
        second.setReadNegativeStrandFlag(true);

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        SAMRecord consensusRead = syncOutcome.CombinedRecord;
        assertNotNull(consensusRead);

        readContextCounter.processRead(consensusRead, 1, new FragmentData(first, second));
        readContextCounter.processRead(consensusRead, 1, new FragmentData(first, second));
        readContextCounter.processRead(consensusRead, 1, new FragmentData(first, second));

        // swap and re-process
        first.setReadNegativeStrandFlag(true);
        second.setReadNegativeStrandFlag(false);

        readContextCounter.processRead(consensusRead, 1, new FragmentData(first, second));

        assertEquals(4, readContextCounter.readStrandBiasAlt().depth());
        assertEquals(0.25, readContextCounter.readStrandBiasAlt().bias(), 0.01);
    }

    @Test
    public void testSplitReads()
    {
        // first a basic match
        SAMRecord first = createSamRecord(
                READ_ID, CHR_1, 20, REF_BASES.substring(20, 30) + REF_BASES.substring(60, 80), "10M30N20M");
        first.setReadNegativeStrandFlag(true);

        SAMRecord second = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 30) + REF_BASES.substring(60, 70), "20M30N10M");

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        SAMRecord combined = syncOutcome.CombinedRecord;
        assertNotNull(combined);
        assertEquals(10, combined.getAlignmentStart());
        assertEquals(79, combined.getAlignmentEnd());
        assertEquals(REF_BASES.substring(10, 30) + REF_BASES.substring(60, 80), combined.getReadString());
        assertEquals("20M30N20M", combined.getCigarString());
    }

    @Test
    public void testMismatches()
    {
        // first a basic match
        SAMRecord first = createSamRecord(
                READ_ID, CHR_1, 1, REF_BASES.substring(1, 11) + "C" + REF_BASES.substring(11, 21), "10M1I10M");

        SAMRecord second = createSamRecord(READ_ID, CHR_1, 5, REF_BASES.substring(5, 25), "20M");
        second.setReadNegativeStrandFlag(true);

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        Assert.assertEquals(FragmentSyncType.CIGAR_MISMATCH, syncOutcome.SyncType);

        // an inversion
        first = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 40), "30M");

        second = createSamRecord(READ_ID, CHR_1, 15, REF_BASES.substring(15, 45), "30M");

        syncOutcome = formFragmentRead(first, second);
        assertEquals(FragmentSyncType.INVERSION, syncOutcome.SyncType);

        // off by 1
        first = createSamRecord(
                READ_ID, CHR_1, 1, REF_BASES.substring(1, 11) + "C" + REF_BASES.substring(11, 21), "10M1I10M");

        second = createSamRecord(
                READ_ID, CHR_1, 2, REF_BASES.substring(1, 12) + "C" + REF_BASES.substring(12, 21), "10M1I10M");
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(FragmentSyncType.CIGAR_MISMATCH, syncOutcome.SyncType);

        /*
        // too many mismatches
        first = createSamRecord(READ_ID, CHR_1, 1, REF_BASES.substring(1, 21), "20M");
        second = createSamRecord(READ_ID, CHR_1, 1, REF_BASES.substring(2, 22), "20M");
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(FragmentSyncType.BASE_MISMATCH, syncOutcome.SyncType);
         */

        // INDEL vs aligned
        first = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 20) + REF_BASES.substring(25, 40), "10M5D15M");

        second = createSamRecord(
                READ_ID, CHR_1, 15, REF_BASES.substring(15, 45), "30M");
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(FragmentSyncType.CIGAR_MISMATCH, syncOutcome.SyncType);

        // different INDELs
        first = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 20) + "C" + REF_BASES.substring(20, 40), "10M1I20M");

        second = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 20) + "CC" + REF_BASES.substring(20, 40), "10M2I20M");
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(FragmentSyncType.CIGAR_MISMATCH, syncOutcome.SyncType);

        // non-overlapping but different INDELs
        first = createSamRecord(
                READ_ID, CHR_1, 1, REF_BASES.substring(1, 6) + "C" + REF_BASES.substring(6, 41), "5M1I35M");

        second = createSamRecord(
                READ_ID, CHR_1, 30, REF_BASES.substring(30, 40) + REF_BASES.substring(45, 75), "10M5D30M");
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(FragmentSyncType.CIGAR_MISMATCH, syncOutcome.SyncType);
    }

    @Test
    public void testTruncatedFragments()
    {
        SAMRecord first = createSamRecord(READ_ID, CHR_1, 5, REF_BASES.substring(5, 35), "30M");

        first.setInferredInsertSize(25);
        first.setReadNegativeStrandFlag(true);

        SAMRecord second = createSamRecord(READ_ID, CHR_1, 10, REF_BASES.substring(10, 40), "30M");
        second.setInferredInsertSize(-25);

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        assertEquals(25, syncOutcome.CombinedRecord.getInferredInsertSize());
        assertEquals(10, syncOutcome.CombinedRecord.getAlignmentStart());
        assertEquals(34, syncOutcome.CombinedRecord.getAlignmentEnd());
        assertEquals("25M", syncOutcome.CombinedRecord.getCigarString());
        assertEquals(25, syncOutcome.CombinedRecord.getBaseQualities().length);

        String readBases = syncOutcome.CombinedRecord.getReadString();
        assertEquals(REF_BASES.substring(10, 35), readBases);
        assertEquals(FragmentSyncType.COMBINED, syncOutcome.SyncType);

        // test again with soft-clips
        // actual alignment is 11 -> 42 = 22 bases
        int fragLength = 22;
        first = createSamRecord(READ_ID, CHR_1, 11, REF_BASES.substring(9, 50), "2S35M4S");
        first.setInferredInsertSize(fragLength);

        second = createSamRecord(READ_ID, CHR_1, 11, REF_BASES.substring(5, 45), "6S32M2S");

        second.setInferredInsertSize(-fragLength);
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(fragLength, syncOutcome.CombinedRecord.getInferredInsertSize());
        assertEquals(11, syncOutcome.CombinedRecord.getAlignmentStart());
        assertEquals(42, syncOutcome.CombinedRecord.getAlignmentEnd());
        assertEquals("2S32M2S", syncOutcome.CombinedRecord.getCigarString());
        assertEquals(36, syncOutcome.CombinedRecord.getBaseQualities().length);

        readBases = syncOutcome.CombinedRecord.getReadString();
        assertEquals(REF_BASES.substring(9, 45), readBases);
        assertEquals(FragmentSyncType.COMBINED, syncOutcome.SyncType);

        // test again with soft-clips on the 3' ends and a DEL

        fragLength = 43 - 8 + 1;
        first = createSamRecord(
                READ_ID, CHR_1, 8, REF_BASES.substring(8, 23) + REF_BASES.substring(28, 45), "15M5D15M2S");
        first.setInferredInsertSize(fragLength);

        second = createSamRecord(
                READ_ID, CHR_1, 11, REF_BASES.substring(5, 23) + REF_BASES.substring(28, 43), "6S12M5D15M");

        second.setInferredInsertSize(-fragLength);
        second.setReadNegativeStrandFlag(true);

        syncOutcome = formFragmentRead(first, second);
        assertEquals(fragLength, syncOutcome.CombinedRecord.getInferredInsertSize());
        assertEquals(8, syncOutcome.CombinedRecord.getAlignmentStart());
        assertEquals(42, syncOutcome.CombinedRecord.getAlignmentEnd());
        assertEquals("15M5D15M", syncOutcome.CombinedRecord.getCigarString());
        assertEquals(30, syncOutcome.CombinedRecord.getBaseQualities().length);

        readBases = syncOutcome.CombinedRecord.getReadString();
        assertEquals(REF_BASES.substring(8, 23) + REF_BASES.substring(28, 43), readBases);
        assertEquals(FragmentSyncType.COMBINED, syncOutcome.SyncType);
    }

    @Test
    public void testTruncatedFragmentIndelMismatch()
    {
        // 3' read has a delete instead of a soft-clip at the start and end
        SAMRecord first = createSamRecord(
                READ_ID, CHR_1, 4, REF_BASES.substring(4, 8) + REF_BASES.substring(10, 40), "4M2D30M");

        first.setInferredInsertSize(30);
        first.setReadNegativeStrandFlag(true);

        SAMRecord second = createSamRecord(
                READ_ID, CHR_1, 10, REF_BASES.substring(10, 40) + REF_BASES.substring(42, 46), "30M2D2M2S");
        second.setInferredInsertSize(-30);

        FragmentSyncOutcome syncOutcome = formFragmentRead(first, second);
        assertEquals(30, syncOutcome.CombinedRecord.getInferredInsertSize());
        assertEquals(10, syncOutcome.CombinedRecord.getAlignmentStart());
        assertEquals(39, syncOutcome.CombinedRecord.getAlignmentEnd());
        assertEquals("30M", syncOutcome.CombinedRecord.getCigarString());
        assertEquals(30, syncOutcome.CombinedRecord.getBaseQualities().length);

        String readBases = syncOutcome.CombinedRecord.getReadString();
        assertEquals(REF_BASES.substring(10, 40), readBases);
        assertEquals(FragmentSyncType.COMBINED, syncOutcome.SyncType);
    }

    private static FragmentSyncOutcome formFragmentRead(final SAMRecord first, final SAMRecord second)
    {
        return CombinedSyncData.formFragmentRead(first, second, null);
    }
}
