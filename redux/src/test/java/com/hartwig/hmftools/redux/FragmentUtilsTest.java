package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.flipFirstInPair;
import static com.hartwig.hmftools.redux.TestUtils.REF_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_CIGAR;
import static com.hartwig.hmftools.redux.old.DuplicateGroupBuilderOld.calcBaseQualAverage;
import static com.hartwig.hmftools.redux.old.DuplicateGroupBuilderOld.findPrimaryFragment;
import static com.hartwig.hmftools.redux.TestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_BASES;
import static com.hartwig.hmftools.redux.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.redux.TestUtils.createFragment;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.redux.common.FragmentCoords;
import com.hartwig.hmftools.redux.old.DuplicateGroupOld;
import com.hartwig.hmftools.redux.old.FragmentOld;
import com.hartwig.hmftools.redux.old.FragmentCoordsOld;
import com.hartwig.hmftools.redux.old.FragmentUtils;
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
        SAMRecord read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, false, null, true, TEST_READ_CIGAR);

        FragmentCoords fragmentCoords = FragmentCoords.fromRead(read);

        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.IsForward);
        assertEquals("1_100_1_299_R", fragmentCoords.Key);

        flipFirstInPair(read);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertFalse(fragmentCoords.IsForward);

        // a mate upper read
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 200, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 100, true, false, null, false, TEST_READ_CIGAR);

        flipFirstInPair(read);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertFalse(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.IsForward);
        assertEquals("1_100_1_299_R", fragmentCoords.Key);

        // reverse orientation fragment with soft-clips at both ends
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, "80M20S",
                CHR_1, 400, true, false, null, false, "10S80M");

        fragmentCoords = FragmentCoords.fromRead(read);
        assertEquals(199, fragmentCoords.PositionLower);
        assertEquals(390, fragmentCoords.PositionUpper);
        assertEquals(REVERSE, fragmentCoords.OrientLower);
        assertEquals(FORWARD, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.IsForward);
        assertEquals("1_199_R_1_390", fragmentCoords.Key);

        // using supplementary data to get primary coords
        SupplementaryReadData suppAlignment = new SupplementaryReadData(
                CHR_1, 100, SupplementaryReadData.SUPP_POS_STRAND, TEST_READ_CIGAR, 60);

        read = createSamRecord(
                TEST_READ_ID, CHR_2, 1000, TEST_READ_BASES, TEST_READ_CIGAR,
                CHR_1, 200, false, true, suppAlignment, true, TEST_READ_CIGAR);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertTrue(fragmentCoords.SupplementarySourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(299, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertEquals(REVERSE, fragmentCoords.OrientUpper);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.IsForward);
        assertEquals("1_100_1_299_R", fragmentCoords.Key);

        // test again with an unmapped mate
        read = createSamRecord(
                TEST_READ_ID, CHR_1, 100, TEST_READ_BASES, TEST_READ_CIGAR,
                NO_CHROMOSOME_NAME, 100, false, false, null);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertFalse(fragmentCoords.UnmappedSourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(NO_POSITION, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.IsForward);
        assertEquals("1_100", fragmentCoords.Key);

        read.setReadNegativeStrandFlag(true);
        fragmentCoords = FragmentCoords.fromRead(read);
        assertEquals(REVERSE, fragmentCoords.OrientLower);
        assertEquals("1_199_R", fragmentCoords.Key);

        // the unmapped read
        read = createSamRecord(
                TEST_READ_ID, NO_CHROMOSOME_NAME, 100, TEST_READ_BASES, NO_CIGAR,
                CHR_1, 100, false, false, null, false, TEST_READ_CIGAR);
        read.setReadUnmappedFlag(true);

        fragmentCoords = FragmentCoords.fromRead(read);
        assertTrue(fragmentCoords.UnmappedSourced);
        assertEquals(100, fragmentCoords.PositionLower);
        assertEquals(NO_POSITION, fragmentCoords.PositionUpper);
        assertEquals(FORWARD, fragmentCoords.OrientLower);
        assertTrue(fragmentCoords.ReadIsLower);
        assertTrue(fragmentCoords.IsForward);
        assertEquals("1_100_U", fragmentCoords.Key);
    }

    @Test
    public void testPrimaryDuplicateIdentification()
    {
        FragmentOld fragment = createFragment(TEST_READ_ID, CHR_1, 100);
        double baseAvg = calcBaseQualAverage(fragment);
        assertEquals(DEFAULT_QUAL, baseAvg, 0.1);

        FragmentOld fragment2 = createFragment(TEST_READ_ID, CHR_1, 100);

        SAMRecord read2 = fragment2.reads().get(0);
        setBaseQualities(read2, DEFAULT_QUAL - 1);

        assertEquals(DEFAULT_QUAL - 1, calcBaseQualAverage(fragment2), 0.1);

        List<FragmentOld> fragments = Lists.newArrayList(fragment2, fragment);

        FragmentOld primary = findPrimaryFragment(fragments, false);
        assertEquals(primary, fragment);
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

        FragmentOld frag1 = new FragmentOld(read1);
        frag1.addRead(mate1);
        assertTrue(frag1.isPreciseInversion());

        SAMRecord read2 = createSamRecord("READ_02", CHR_1, posStart, TEST_READ_BASES, "100M", CHR_1, posStart,
                false, false, null);
        read2.setFirstOfPairFlag(true);

        SAMRecord mate2 = createSamRecord("READ_02", CHR_1, posStart, TEST_READ_BASES, "100M", CHR_1, posStart,
                false, false, null);
        mate2.setSecondOfPairFlag(true);
        mate2.setFirstOfPairFlag(false);

        FragmentOld frag2 = new FragmentOld(read2);
        frag2.addRead(mate2);
        assertTrue(frag2.isPreciseInversion());

        DuplicateGroupOld duplicateGroup = new DuplicateGroupOld(null, frag1);
        duplicateGroup.addFragment(frag2);

        duplicateGroup.categoriseReads();

        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES);
        ConsensusReads consensusReads = new ConsensusReads(refGenome);
        List<SAMRecord> reads = duplicateGroup.popCompletedReads(consensusReads, false);
        assertEquals(6, reads.size());
    }
}
