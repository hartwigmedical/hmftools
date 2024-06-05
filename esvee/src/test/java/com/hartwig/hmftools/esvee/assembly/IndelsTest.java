package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.region.Orientation;
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
