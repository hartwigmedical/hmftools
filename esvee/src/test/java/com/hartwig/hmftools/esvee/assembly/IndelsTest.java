package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_RANDOM_100;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.calcIndelInferredUnclippedPositions;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import org.junit.Test;

import htsjdk.samtools.CigarElement;

public class IndelsTest
{
    @Test
    public void testIndelCoords()
    {
        int readStart = 100;

        String cigar = "20M2I20M";
        List<CigarElement> cigarElements = CigarUtils.cigarElementsFromStr(cigar);

        // too short
        IndelCoords indelCoords = findIndelCoords(readStart, cigarElements, 10);
        assertNull(indelCoords);

        // valid delete
        cigar = "20M10D20M"; // 100-119 M, 120-129 D, 130-139 M
        cigarElements = CigarUtils.cigarElementsFromStr(cigar);
        indelCoords = findIndelCoords(readStart, cigarElements, 10);
        assertNotNull(indelCoords);
        assertEquals(10, indelCoords.Length);
        assertEquals(119, indelCoords.PosStart);
        assertEquals(130, indelCoords.PosEnd);

        // valid insert
        cigar = "20M10I20M"; // 100-119 M, 10 I @ 119, 120-129 M
        cigarElements = CigarUtils.cigarElementsFromStr(cigar);
        indelCoords = findIndelCoords(readStart, cigarElements, 10);
        assertNotNull(indelCoords);
        assertEquals(10, indelCoords.Length);
        assertEquals(119, indelCoords.PosStart);
        assertEquals(120, indelCoords.PosEnd);

        // multiple, chooses the longest
        cigar = "10M10I10M15D10M20I10M25D20M"; // 100-109 M, 10 I @ 109 110-119 M, 120-134 D, 135-144M, 20 I @ 144, 145-154M, 155-179 D, 180-199M
        cigarElements = CigarUtils.cigarElementsFromStr(cigar);
        indelCoords = findIndelCoords(readStart, cigarElements, 10);
        assertNotNull(indelCoords);
        assertEquals(25, indelCoords.Length);
        assertEquals(154, indelCoords.PosStart);
        assertEquals(180, indelCoords.PosEnd);
    }

    @Test
    public void testImpliedEdgeIndelsToSoftClip()
    {
        // the cigar remains unch but the implied new alignments are calculated and the unclipped alignments stored

        // converts both sides of the insert
        String cigar = "20M10I20M";
        String readBases = REF_BASES_RANDOM_100.substring(0, 50);
        Read read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertTrue(calcIndelInferredUnclippedPositions(read));
        assertEquals(100, read.alignmentStart());
        assertEquals(139, read.alignmentEnd());
        assertEquals(100, read.unclippedStart());
        assertEquals(139, read.unclippedEnd());

        assertTrue(read.hasIndelImpliedUnclippedEnd());
        assertTrue(read.hasIndelImpliedUnclippedStart());
        assertEquals(149, read.indelImpliedUnclippedEnd());
        assertEquals(90, read.indelImpliedUnclippedStart());

        // no conversion for soft-clipped reads
        cigar = "5S20M10I20M5S";
        readBases = REF_BASES_RANDOM_100.substring(0, 60);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(read.hasIndelImpliedUnclippedEnd());
        assertFalse(read.hasIndelImpliedUnclippedStart());

        // no conversion for short indels
        cigar = "20M2I20M2D20M";
        readBases = REF_BASES_RANDOM_100.substring(0, 64);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(read.hasIndelImpliedUnclippedEnd());
        assertFalse(read.hasIndelImpliedUnclippedStart());

        // converts both sides of the delete
        cigar = "10M20D20M10D10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 40);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(169, read.alignmentEnd());

        assertTrue(calcIndelInferredUnclippedPositions(read));
        assertEquals(139, read.indelImpliedUnclippedEnd());
        assertEquals(130, read.indelImpliedUnclippedStart());
    }

    /*
    private static boolean convertEdgeIndelsToSoftClip(final Read read, boolean allowDoubleConversion)
    {
        return ReadAdjustments.convertEdgeIndelsToSoftClip(read, 6, 15);
    }

    @Test
    public void testImpliedEdgeIndelsToSoftClip()
    {
        // the cigar remains unch but the implied new alignments are calculated and the unclipped alignments stored

        // converts both sides of the insert
        String cigar = "18M6I17M";
        String readBases = REF_BASES_RANDOM_100.substring(0, 40);
        Read read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertTrue(convertEdgeIndelsToSoftClip(read, true));
        assertEquals(100, read.alignmentStart());
        assertEquals(134, read.alignmentEnd());
        assertEquals(100, read.unclippedStart());
        assertEquals(134, read.unclippedEnd());
        assertEquals(118, read.indelImpliedAlignmentStart());
        assertEquals(117, read.indelImpliedAlignmentEnd());
        assertEquals(94, read.minUnclippedStart());
        assertEquals(140, read.maxUnclippedEnd());

        // indel too short
        cigar = "10M5I10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 25);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertFalse(convertEdgeIndelsToSoftClip(read, true));

        // converts both sides of the delete
        cigar = "10M6D15M";
        readBases = REF_BASES_RANDOM_100.substring(0, 31);
        read = createRead(TEST_READ_ID, 100, readBases, cigar);
        assertEquals(130, read.alignmentEnd());
        assertTrue(convertEdgeIndelsToSoftClip(read, true));
        assertEquals(100, read.alignmentStart());
        assertEquals(130, read.alignmentEnd());
        assertEquals(100, read.unclippedStart());
        assertEquals(130, read.unclippedEnd());
        assertEquals(116, read.indelImpliedAlignmentStart());
        assertEquals(109, read.indelImpliedAlignmentEnd());
        assertEquals(100, read.minUnclippedStart());
        assertEquals(130, read.maxUnclippedEnd());

        cigar = "70M20D10M";
        readBases = REF_BASES_RANDOM_100.substring(0, 80);
        read = createRead(TEST_READ_ID, 101, readBases, cigar);
        assertEquals(200, read.alignmentEnd());
        assertEquals(200, read.unclippedEnd());
        assertTrue(ReadAdjustments.convertEdgeIndelsToSoftClip(read));
        // assertEquals(100, read.alignmentStart());
        //assertEquals(100, read.unclippedStart());

        // unch
        assertEquals(200, read.alignmentEnd());
        assertEquals(200, read.unclippedEnd());

        // assertEquals(116, read.indelImpliedAlignmentStart());
        assertEquals(170, read.indelImpliedAlignmentEnd());
        //assertEquals(100, read.minUnclippedStart());
        assertEquals(200, read.maxUnclippedEnd());
    }
    */

    @Test
    public void testIndelReadsSupportingSplitJunctions()
    {
        String extBases = REF_BASES_400.substring(300, 350);

        Junction posJunction = new Junction(CHR_1, 100, FORWARD);

        // assembly has 2 high-qual junction reads
        String readBases = REF_BASES_400.substring(51, 101) + extBases;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases, makeCigarString(readBases, 0, extBases.length()));

        String readBases2 = readBases.substring(1);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 52, readBases, makeCigarString(readBases2, 0, extBases.length()));

        // various indel reads that support the junction
        Read indel1 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases, "40M10I50M");
        calcIndelInferredUnclippedPositions(indel1);

        Read indel2 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases, "40M10D60M");
        calcIndelInferredUnclippedPositions(indel2);

        List<Read> reads = List.of(read1, read2, indel1, indel2);

        JunctionAssembler junctionAssembler = new JunctionAssembler(posJunction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(4, assembly.supportCount());
        assertEquals(0, assembly.mismatchReadCount());

        // note the indel with 10D does not support the ref but this not currently checked
    }

    @Test
    public void testLongDeleteAssemblies()
    {
        Junction posJunction = new Junction(CHR_1, 50, FORWARD, false, true, false);

        // first a basic assembly with all reads agreeing
        String readBases = REF_BASES_200.substring(11, 51) + REF_BASES_200.substring(100, 140);
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "40M49D40M");

        readBases = REF_BASES_200.substring(21, 51) + REF_BASES_200.substring(100, 150);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 21, readBases, "30M49D50M");

        // other reads will soft-clip at the junctions
        readBases = REF_BASES_200.substring(11, 51) + REF_BASES_200.substring(100, 120);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "40M20S");

        List<Read> reads = List.of(read1, read2, read3);

        JunctionAssembler junctionAssembler = new JunctionAssembler(posJunction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(3, assembly.supportCount());
        assertEquals(0, assembly.mismatchReadCount());

        // test the other side
        Junction negJunction = new Junction(CHR_1, 100, REVERSE, false, true, false);

        readBases = REF_BASES_200.substring(31, 51) + REF_BASES_200.substring(100, 160);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, "20S60M");

        reads = List.of(read1, read2, read3);

        junctionAssembler = new JunctionAssembler(negJunction);
        assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        assembly = assemblies.get(0);
        assertEquals(3, assembly.supportCount());
    }
}
