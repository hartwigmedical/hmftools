package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.calcIndelInferredUnclippedPositions;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readSoftClipsAndCrossesJunction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

import org.junit.Test;

public class ReadUtilsTest
{
    @Test
    public void testForwardReadParseState()
    {
        String readBases = REF_BASES_200.substring(1, 41);
        Read read = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "10S30M");
        int readJunctionIndex = 10; // soft-clip length
        boolean moveForward = true;

        ReadParseState readState = new ReadParseState(moveForward, read, readJunctionIndex, true);

        assertEquals(11, readState.refPosition());
        assertEquals(10, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(11), (char)readState.currentBase());

        for(int i = 0; i < 29; ++i)
        {
            readState.moveNext();
        }

        assertEquals(40, readState.refPosition());
        assertEquals(39, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(40), (char)readState.currentBase());
        assertFalse(readState.exhausted());

        readState.moveNext();
        assertTrue(readState.exhausted());

        // test again with indels
        read = createRead(READ_ID_GENERATOR.nextId(), 1, REF_BASES_200.substring(1, 36), "5M5D5M5I10M10S");

        readState = new ReadParseState(moveForward, read, 0, true);

        assertEquals(1, readState.refPosition());
        assertEquals(0, readState.readIndex());
        assertEquals(M, readState.operator());

        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(6, readState.refPosition());
        assertEquals(4, readState.readIndex());
        assertEquals(D, readState.operator());

        // move through the delete
        readState.moveNext();

        assertEquals(7, readState.refPosition());
        assertEquals(4, readState.readIndex());
        assertEquals(D, readState.operator());

        readState.moveNext();
        readState.moveNext();
        readState.moveNext();

        assertEquals(10, readState.refPosition());
        assertEquals(4, readState.readIndex());
        assertEquals(D, readState.operator());

        readState.moveNext();

        assertEquals(11, readState.refPosition());
        assertEquals(5, readState.readIndex());
        assertEquals(M, readState.operator());

        // move onto insert
        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(15, readState.refPosition());
        assertEquals(10, readState.readIndex());
        assertEquals(I, readState.operator());

        readState.moveNext();
        assertEquals(15, readState.refPosition());
        assertEquals(11, readState.readIndex());
        assertEquals(I, readState.operator());

        // through rest of insert
        readState.moveNext();
        readState.moveNext();
        readState.moveNext();
        readState.moveNext();

        assertEquals(16, readState.refPosition());
        assertEquals(15, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(16), (char)readState.currentBase());

        // stop at final soft-clipping
        for(int i = 0; i < 10; ++i)
        {
            readState.moveNext();
        }

        assertTrue(readState.exhausted());
    }

    @Test
    public void testReverseReadState()
    {
        String readBases = REF_BASES_200.substring(1, 41);
        Read read = createRead(READ_ID_GENERATOR.nextId(), 1, readBases, "30M10S");
        int readJunctionIndex = 29;
        boolean moveForward = false;

        ReadParseState readState = new ReadParseState(moveForward, read, readJunctionIndex, true);

        assertEquals(30, readState.refPosition());
        assertEquals(29, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(30), (char)readState.currentBase());

        for(int i = 0; i < 29; ++i)
        {
            readState.moveNext();
        }

        assertEquals(1, readState.refPosition());
        assertEquals(0, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(1), (char)readState.currentBase());
        assertFalse(readState.exhausted());

        readState.moveNext();
        assertTrue(readState.exhausted());

        // test again with indels
        read = createRead(READ_ID_GENERATOR.nextId(), 11, REF_BASES_200.substring(1, 36), "10S10M5D5M5I5M");
        // bases 31-35, index 30-34 = 5M
        // bases 30, index 25-29 = 5I
        // bases 26-30, index 20-24 = 5M
        // bases 21-25,  = 5D
        // bases 11-20, index 10-19 = 10M
        // bases 1-10, index 0-9 = 10S

        readJunctionIndex = read.basesLength() - 1;
        readState = new ReadParseState(moveForward, read, readJunctionIndex, true);

        assertEquals(35, readState.refPosition());
        assertEquals(34, readState.readIndex());
        assertEquals(M, readState.operator());

        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(31, readState.refPosition());
        assertEquals(29, readState.readIndex());
        assertEquals(I, readState.operator());

        // move through the insert
        readState.moveNext();

        assertEquals(31, readState.refPosition());
        assertEquals(28, readState.readIndex());
        assertEquals(I, readState.operator());

        readState.moveNext();
        readState.moveNext();
        readState.moveNext();
        readState.moveNext();

        assertEquals(30, readState.refPosition());
        assertEquals(24, readState.readIndex());
        assertEquals(M, readState.operator());

        // move to start of delete
        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(25, readState.refPosition());
        assertEquals(20, readState.readIndex());
        assertEquals(D, readState.operator());

        // move past delete
        for(int i = 0; i < 5; ++i)
        {
            readState.moveNext();
        }

        assertEquals(20, readState.refPosition());
        assertEquals(19, readState.readIndex());
        assertEquals(M, readState.operator());

        // move to last aligned base
        for(int i = 0; i < 9; ++i)
        {
            readState.moveNext();
        }

        assertEquals(11, readState.refPosition());
        assertEquals(10, readState.readIndex());
        assertEquals(M, readState.operator());

        // stop at final soft-clipping
        readState.moveNext();
        assertTrue(readState.exhausted());
    }

    @Test
    public void testHardClippedReads()
    {
        // index 0-14, ref pos 1-15
        // index DEL, ref pos 16-21
        // index 15-46, ref pos 22-52
        String cigar = "105H15M6D31M";
        String readBases = REF_BASES_200.substring(1, 16) + REF_BASES_200.substring(22, 53);

        Read read = createRead(READ_ID_GENERATOR.nextId(), 1, readBases, cigar);
        int readJunctionIndex = 15;
        boolean moveForwards = true;

        ReadParseState readState = new ReadParseState(moveForwards, read, readJunctionIndex, true);

        assertEquals(22, readState.refPosition());
        assertEquals(15, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(22), (char)readState.currentBase());

        readState.moveNext();

        assertEquals(23, readState.refPosition());
        assertEquals(16, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(23), (char)readState.currentBase());

        // index 0-30, ref pos 1-31
        // index DEL, ref pos 32-37
        // index 31-45, ref pos 38-52
        cigar = "31M6D15M105H";
        readBases = REF_BASES_200.substring(1, 32) + REF_BASES_200.substring(38, 53);

        read = createRead(READ_ID_GENERATOR.nextId(), 1, readBases, cigar);
        readJunctionIndex = 31;
        moveForwards = false;

        readState = new ReadParseState(moveForwards, read, readJunctionIndex, true);

        assertEquals(38, readState.refPosition());
        assertEquals(31, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(38), (char)readState.currentBase());

        readState.moveNext();

        assertEquals(37, readState.refPosition());
        assertEquals(31, readState.readIndex());
        assertEquals(D, readState.operator());

        for(int i = 0; i < 6; ++i)
        {
            readState.moveNext();
        }

        assertEquals(31, readState.refPosition());
        assertEquals(30, readState.readIndex());
        assertEquals(M, readState.operator());
        assertEquals(REF_BASES_200.charAt(31), (char)readState.currentBase());
    }

    @Test
    public void testMiscReadFunctions()
    {
        String cigarStr = "5S20M8S";
        Read read = createRead(TEST_READ_ID, 20, REF_BASES_RANDOM_100.substring(15, 48), cigarStr);

        assertEquals(15, read.unclippedStart());
        assertEquals(47, read.unclippedEnd());
        assertEquals(3, read.cigarElements().size());
        assertEquals(cigarStr, read.cigarString());

        assertEquals(5, read.getReadIndexAtReferencePosition(20));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(19));
        assertEquals(4, read.getReadIndexAtReferencePosition(19, true));
        assertEquals(0, read.getReadIndexAtReferencePosition(15, true));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(14, true));

        assertEquals(24, read.getReadIndexAtReferencePosition(39));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(40));
        assertEquals(25, read.getReadIndexAtReferencePosition(40, true));
        assertEquals(32, read.getReadIndexAtReferencePosition(47, true));
        assertEquals(INVALID_INDEX, read.getReadIndexAtReferencePosition(48, true));

        // with INDELs

        cigarStr = "4S20M10D10M5I10M3S";
        // positions 4S then 20-39 = M,   40-49 = D, 50-59 = M,    T at 59,     60-69 = M, then 3S
        // index: 4S is 0-3, 20M is 4-23, D is none, 10M is 24-33, 5I is 34-38, 10M is 39-48, 3S is 49-51
        read = createRead(
                TEST_READ_ID, 20,
                REF_BASES_RANDOM_100.substring(16, 40) + REF_BASES_RANDOM_100.substring(50, 60) + "GGGGG" + REF_BASES_RANDOM_100.substring(60, 73),
                cigarStr);

        assertEquals(16, read.unclippedStart());
        assertEquals(69, read.alignmentEnd());
        assertEquals(72, read.unclippedEnd());
        assertEquals(7, read.cigarElements().size());
        assertEquals(52, read.basesLength());

        assertEquals(0, read.getReadIndexAtReferencePosition(16, true));
        assertEquals(4, read.getReadIndexAtReferencePosition(20, false));
        assertEquals(5, read.getReadIndexAtReferencePosition(21, false));
        assertEquals(23, read.getReadIndexAtReferencePosition(39, false)); // start of DEL
        assertEquals(23, read.getReadIndexAtReferencePosition(48, false)); // within DEL

        assertEquals(24, read.getReadIndexAtReferencePosition(50, false)); // end of DEL
        assertEquals(33, read.getReadIndexAtReferencePosition(59, false)); // end of next M
        assertEquals(39, read.getReadIndexAtReferencePosition(60, false)); // base after insert
        assertEquals(48, read.getReadIndexAtReferencePosition(69, false)); // final aligned base
        assertEquals(49, read.getReadIndexAtReferencePosition(70, true));
        assertEquals(50, read.getReadIndexAtReferencePosition(71, true));
        assertEquals(51, read.getReadIndexAtReferencePosition(72, true)); // final SC base

        // compare to SAM record method, which is 1-based
        assertEquals(5, read.bamRecord().getReadPositionAtReferencePosition(20));
        assertEquals(24, read.bamRecord().getReadPositionAtReferencePosition(39)); // start of DEL

        assertEquals(25, read.bamRecord().getReadPositionAtReferencePosition(50)); // end of DEL
        assertEquals(34, read.bamRecord().getReadPositionAtReferencePosition(59)); // end of next M
        assertEquals(40, read.bamRecord().getReadPositionAtReferencePosition(60)); // base after insert
        assertEquals(49, read.bamRecord().getReadPositionAtReferencePosition(69)); // final aligned base
    }

    @Test
    public void testReadJunctionConditions()
    {
        Junction posJunction = new Junction(CHR_1, 30, FORWARD);
        String readBases = REF_BASES_200.substring(10, 20);
        String readId = "READ_01";

        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_200.substring(1, 100));

        Read read = createRead(readId, 21, readBases, "10M10S");
        assertTrue(readSoftClipsAndCrossesJunction(read, posJunction, refGenome));

        read = createRead(readId, 23, readBases, "10M10S");
        assertTrue(readSoftClipsAndCrossesJunction(read, posJunction, refGenome));

        read = createRead(readId, 19, readBases, "10M10S");
        assertTrue(readSoftClipsAndCrossesJunction(read, posJunction, refGenome));

        // alignment too far from junction
        read = createRead(readId, 24, readBases, "10M10S");
        assertFalse(readSoftClipsAndCrossesJunction(read, posJunction, refGenome));

        // any realigned indel that crosses the junction is permitted
        read = createRead(readId, 10, readBases, "10M10I10M");
        assertTrue(calcIndelInferredUnclippedPositions(read));
        assertTrue(readSoftClipsAndCrossesJunction(read, posJunction, refGenome));

        // a read with a low count of SNVs
        readBases = REF_BASES_200.substring(20, 30) + "GGGGG";
        read = createRead(readId, 20, readBases, "15M");
        read.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(readSoftClipsAndCrossesJunction(read, posJunction, refGenome));

        Junction negJunction = new Junction(CHR_1, 30, REVERSE);

        read = createRead(readId, 30, readBases, "10S10M");
        assertTrue(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));

        read = createRead(readId, 28, readBases, "10S10M");
        assertTrue(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));

        read = createRead(readId, 32, readBases, "10S10M");
        assertTrue(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));

        // alignment too far from junction
        read = createRead(readId, 33, readBases, "10S10M");
        assertFalse(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));

        // any realigned indel that crosses the junction is permitted
        read = createRead(readId, 10, readBases, "10M10I10M");
        assertTrue(calcIndelInferredUnclippedPositions(read));
        assertTrue(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));


        // a read with a low count of SNVs
        readBases = "GGGG" + REF_BASES_200.substring(30, 50);
        read = createRead(readId, 26, readBases, "34M");
        read.bamRecord().setAttribute(NUM_MUTATONS_ATTRIBUTE, 2);
        assertTrue(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));

        readBases = REF_BASES_200.substring(27, 51);
        read = createRead(readId, 26, readBases, "34M");
        assertFalse(readSoftClipsAndCrossesJunction(read, negJunction, refGenome));
    }

}
