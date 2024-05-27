package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.common.DuplicateGroupBuilder.calcBaseQualAverage;
import static com.hartwig.hmftools.redux.common.DuplicateGroupBuilder.findPrimaryFragment;
import static com.hartwig.hmftools.redux.common.FragmentCoordinates.NO_COORDS;
import static com.hartwig.hmftools.redux.common.FragmentUtils.formChromosomePartition;
import static com.hartwig.hmftools.redux.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.TestUtils.createFragment;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.common.DuplicateGroup;
import com.hartwig.hmftools.redux.common.Fragment;
import com.hartwig.hmftools.redux.common.FragmentCoordinates;
import com.hartwig.hmftools.redux.common.FragmentUtils;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import static com.hartwig.hmftools.redux.TestUtils.setBaseQualities;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class FragmentUtilsTest
{
    @Test
    public void testFragmentCoords()
    {
        SAMRecord read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                false, false, null);

        int ucPos = getFivePrimeUnclippedPosition(read);
        assertEquals(100, ucPos);

        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "5S95M", CHR_1, 200,
                false, false, null);

        ucPos = getFivePrimeUnclippedPosition(read);
        assertEquals(95, ucPos);

        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "5S80M15S", CHR_1, 200,
                true, false, null);

        ucPos = getFivePrimeUnclippedPosition(read);
        assertEquals(194, ucPos);

        // test a read pair with one read unmapped

        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "*", CHR_1, 100,
                false, false, null);

        read.setReadUnmappedFlag(true);
        read.setInferredInsertSize(100);

        Fragment fragment = new Fragment(read);

        assertEquals(NO_COORDS, fragment.coordinates());

        SAMRecord mateRead = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR, CHR_1, 100,
                false, false, null);

        mateRead.setMateUnmappedFlag(true);
        mateRead.setInferredInsertSize(-100);

        fragment.addRead(mateRead);

        assertEquals("1_100_100", fragment.coordinates().Key);
        assertFalse(fragment.coordinates().Incomplete);
    }

    private static FragmentCoordinates getFragmentCoordinates(final SAMRecord read)
    {
        return FragmentUtils.getFragmentCoordinates(Lists.newArrayList(read), true);
    }

    @Test
    public void testFragmentCoordsFromCigars()
    {
        // test coordinates from CIGAR strings
        String cigarStr = "5S100M5S";
        int ucPos = getFivePrimeUnclippedPosition(100, cigarStr, true);
        assertEquals(95, ucPos);

        ucPos = getFivePrimeUnclippedPosition(100, cigarStr, false);
        assertEquals(204, ucPos);

        cigarStr = "5S10M5D10I15D20M5S"; // +10 + 5, ignore 10, +15, +20, +5 = 40
        ucPos = getFivePrimeUnclippedPosition(100, cigarStr, false);
        assertEquals(154, ucPos);

        // RNA
        cigarStr = "5S60M2000N40M10S";
        ucPos = getFivePrimeUnclippedPosition(100, cigarStr, true); // -5
        assertEquals(95, ucPos);

        ucPos = getFivePrimeUnclippedPosition(100, cigarStr, false); // +60, +2000, +40, +10
        assertEquals(2209, ucPos);

        // test coordinates from reads
        SAMRecord read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                false, false, null);
        read.setMateNegativeStrandFlag(true);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        FragmentCoordinates fragmentCoords = getFragmentCoordinates(read);
        assertEquals("1_100_1_299_R", fragmentCoords.Key);
        assertEquals(100, fragmentCoords.InitialPosition);

        // mate on earlier chromosome
        read = createSamRecord(TEST_READ_ID, CHR_2, 100, TEST_READ_BASES, "100M", CHR_1, 200,
                false, false, null);
        read.setMateNegativeStrandFlag(false);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        fragmentCoords = getFragmentCoordinates(read);
        assertEquals("1_200_2_100_N", fragmentCoords.keyOriented());
        assertEquals(200, fragmentCoords.InitialPosition);

        // mate in earlier position, and fragment reversed
        read = createSamRecord(TEST_READ_ID, CHR_1, 200, TEST_READ_BASES, "100M", CHR_1, 100,
                false, false, null);
        read.setMateNegativeStrandFlag(true);
        read.setSecondOfPairFlag(true);
        read.setFirstOfPairFlag(false);
        read.setAttribute(MATE_CIGAR_ATTRIBUTE, "100M");

        fragmentCoords = getFragmentCoordinates(read);
        assertEquals("1_199_R_1_200", fragmentCoords.Key);
        assertEquals(-199, fragmentCoords.InitialPosition);

        // unmapped mate
        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "100M", "", 0,
                false, false, null);
        read.setInferredInsertSize(400);

        fragmentCoords = getFragmentCoordinates(read);
        assertEquals("1_100_400", fragmentCoords.Key);
        assertEquals(100, fragmentCoords.InitialPosition);

        // fragment reversed
        read = createSamRecord(TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "100M", "", 0,
                true, false, null);
        read.setInferredInsertSize(400);

        fragmentCoords = getFragmentCoordinates(read);
        assertEquals("1_199_R_400", fragmentCoords.Key);
        assertEquals(-199, fragmentCoords.InitialPosition);

        // missing mate CIGAR attribute
        read = createSamRecord(TEST_READ_ID, CHR_1, 200, TEST_READ_BASES, "100M", CHR_1, 100,
                false, false, null);

        fragmentCoords = getFragmentCoordinates(read);
        assertEquals("1_200", fragmentCoords.Key);
        assertTrue(fragmentCoords.Incomplete);
    }

    @Test
    public void testPrimaryDuplicateIdentification()
    {
        Fragment fragment = createFragment(TEST_READ_ID, CHR_1, 100);
        double baseAvg = calcBaseQualAverage(fragment);
        assertEquals(DEFAULT_QUAL, baseAvg, 0.1);

        Fragment fragment2 = createFragment(TEST_READ_ID, CHR_1, 100);

        SAMRecord read2 = fragment2.reads().get(0);
        setBaseQualities(read2, DEFAULT_QUAL - 1);

        assertEquals(DEFAULT_QUAL - 1, calcBaseQualAverage(fragment2), 0.1);

        List<Fragment> fragments = Lists.newArrayList(fragment2, fragment);

        Fragment primary = findPrimaryFragment(fragments, false);
        assertEquals(primary, fragment);
    }

    @Test
    public void testReadChromosomePartition()
    {
        assertEquals("1_0", formChromosomePartition(CHR_1, 1, 1000));
        assertEquals("1_0", formChromosomePartition(CHR_1, 999, 1000));
        assertEquals("1_1", formChromosomePartition(CHR_1, 1000, 1000));
    }

    @Test
    public void testDuplicateGroupReadAllocation()
    {
        int posStart = 1;
        SAMRecord read1 = createSamRecord("READ_01", CHR_1, posStart, TEST_READ_BASES, "100M", CHR_1, posStart,
                false, false, null);
        read1.setFirstOfPairFlag(true);

        SAMRecord mate1 = createSamRecord("READ_01", CHR_1, posStart, TEST_READ_BASES, "100M", CHR_1, posStart,
                false, false, null);
        mate1.setSecondOfPairFlag(true);
        mate1.setFirstOfPairFlag(false);

        Fragment frag1 = new Fragment(read1);
        frag1.addRead(mate1);
        assertTrue(frag1.isPreciseInversion());

        SAMRecord read2 = createSamRecord("READ_02", CHR_1, posStart, TEST_READ_BASES, "100M", CHR_1, posStart,
                false, false, null);
        read2.setFirstOfPairFlag(true);

        SAMRecord mate2 = createSamRecord("READ_02", CHR_1, posStart, TEST_READ_BASES, "100M", CHR_1, posStart,
                false, false, null);
        mate2.setSecondOfPairFlag(true);
        mate2.setFirstOfPairFlag(false);

        Fragment frag2 = new Fragment(read2);
        frag2.addRead(mate2);
        assertTrue(frag2.isPreciseInversion());

        DuplicateGroup duplicateGroup = new DuplicateGroup(null, frag1);
        duplicateGroup.addFragment(frag2);

        duplicateGroup.categoriseReads();

        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        ConsensusReads consensusReads = new ConsensusReads(refGenome);
        List<SAMRecord> reads = duplicateGroup.popCompletedReads(consensusReads, false);
        assertEquals(6, reads.size());

    }
}
