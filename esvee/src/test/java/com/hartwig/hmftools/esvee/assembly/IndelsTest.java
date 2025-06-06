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
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO;
import static com.hartwig.hmftools.esvee.assembly.AssemblyDeduper.dedupProximateAssemblies;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.calcIndelInferredUnclippedPositions;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarUtils;
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

        // other reads will soft-clip at or near the junctions
        readBases = REF_BASES_200.substring(11, 51) + REF_BASES_200.substring(100, 120);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "40M20S");

        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "45M15S");
        Read read5 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "35M25S");

        List<Read> reads = List.of(read1, read2, read3, read4, read5);

        JunctionAssembler junctionAssembler = new JunctionAssembler(posJunction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(5, assembly.supportCount());
        assertEquals(0, assembly.mismatchReadCount());

        // test the other side
        Junction negJunction = new Junction(CHR_1, 100, REVERSE, false, true, false);

        readBases = REF_BASES_200.substring(31, 51) + REF_BASES_200.substring(100, 160);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, "20S60M");
        read4 = createRead(READ_ID_GENERATOR.nextId(), 105, readBases, "25S55M");

        reads = List.of(read1, read2, read3, read4);

        junctionAssembler = new JunctionAssembler(negJunction);
        assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        assembly = assemblies.get(0);
        assertEquals(4, assembly.supportCount());
    }

    @Test
    public void testIndelAssemblyDedup()
    {
        Junction indelJunctionPos1 = new Junction(CHR_1, 130, FORWARD, false, true, false);
        Junction indelJunctionNeg1 = new Junction(CHR_1, 131, REVERSE, false, true, false);

        String leftRefBases = REF_BASES_400.substring(100, 130);
        String rightRefBases = REF_BASES_400.substring(131, 171);
        String insertBases = REF_BASES_400.substring(200, 240);

        String indelReadBases = leftRefBases + insertBases + rightRefBases;

        // both pairs of junctions share a read which supports them both, so are not deduped
        Read indelRead = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 101,
                indelReadBases,
                "30M40I30M", CHR_1, 300, true);

        JunctionAssembly indelAssemblyPos1 = new JunctionAssembly(
                indelJunctionPos1, indelRead.getBases(), indelRead.getBaseQuality(), leftRefBases.length());

        indelAssemblyPos1.addJunctionRead(indelRead);
        indelAssemblyPos1.setIndelCoords(indelRead.indelCoords());

        JunctionAssembly indelAssemblyNeg1 = new JunctionAssembly(
                indelJunctionNeg1, indelRead.getBases(), indelRead.getBaseQuality(),
                leftRefBases.length() + insertBases.length() + 1);

        indelAssemblyNeg1.addJunctionRead(indelRead);
        indelAssemblyNeg1.setIndelCoords(indelRead.indelCoords());

        // create a pair where the pos assemblies don't match and the negatives do
        Junction indelJunctionPos2 = new Junction(CHR_1, 140, FORWARD, false, true, false);
        Junction indelJunctionNeg2 = new Junction(CHR_1, 141, REVERSE, false, true, false);

        String insertBases2 = REF_BASES_400.substring(210, 255);

        String indelReadBases2 = leftRefBases + insertBases2 + rightRefBases;

        Read indelRead2Pos = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 111,
                indelReadBases2,
                "30M45I30M", CHR_1, 300, true);

        JunctionAssembly indelAssemblyPos2 = new JunctionAssembly(
                indelJunctionPos2, indelRead2Pos.getBases(), indelRead2Pos.getBaseQuality(), leftRefBases.length());

        indelAssemblyPos2.addJunctionRead(indelRead2Pos);
        indelAssemblyPos2.setIndelCoords(indelRead2Pos.indelCoords());

        Read indelRead2Neg = createRead(
                READ_ID_GENERATOR.nextId(), CHR_1, 111,
                indelReadBases2,
                "30M45I30M", CHR_1, 300, true);

        JunctionAssembly indelAssemblyNeg2 = new JunctionAssembly(
                indelJunctionNeg2, indelRead2Neg.getBases(), indelRead2Neg.getBaseQuality(),
                leftRefBases.length() + insertBases.length() + 1);

        indelAssemblyNeg2.addJunctionRead(indelRead2Neg);
        indelAssemblyNeg2.setIndelCoords(indelRead2Neg.indelCoords());

        // test 1 - pos assemblies don't match, neg assemblies do - first (existing) neg assembly is chosen so both second assemblies are dropped
        List<JunctionAssembly> existingAssemblies = Lists.newArrayList(indelAssemblyPos1);
        List<JunctionAssembly> newAssemblies = Lists.newArrayList(indelAssemblyPos2);
        List<JunctionAssembly> dedupIndels = Lists.newArrayList();

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(1, existingAssemblies.size());
        assertEquals(1, newAssemblies.size());

        existingAssemblies.add(indelAssemblyPos2);

        replicateSupport(indelAssemblyNeg1, indelRead, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO); // ensures this assembly is kept
        newAssemblies = Lists.newArrayList(indelAssemblyNeg1);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(2, existingAssemblies.size());
        assertEquals(1, newAssemblies.size());

        existingAssemblies.add(indelAssemblyNeg1);

        newAssemblies = Lists.newArrayList(indelAssemblyNeg2);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(2, existingAssemblies.size());
        assertFalse(existingAssemblies.contains(indelAssemblyPos2));
        assertEquals(0, newAssemblies.size());
        assertTrue(dedupIndels.isEmpty());

        // test 2 - last neg assembly is kept, first pair is removed
        indelAssemblyNeg1.support().clear();
        indelAssemblyNeg1.addJunctionRead(indelRead);

        replicateSupport(indelAssemblyNeg2, indelRead2Neg, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO); // ensures this assembly is kept

        existingAssemblies = Lists.newArrayList(indelAssemblyPos1, indelAssemblyPos2, indelAssemblyNeg1);

        newAssemblies = Lists.newArrayList(indelAssemblyNeg2);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(1, existingAssemblies.size());
        assertTrue(existingAssemblies.contains(indelAssemblyPos2));
        assertFalse(existingAssemblies.contains(indelAssemblyPos1));
        assertFalse(existingAssemblies.contains(indelAssemblyNeg1));
        assertEquals(1, newAssemblies.size());
        assertTrue(dedupIndels.isEmpty());

        // test 3 - positives are deduped first, neg deduped based on this
        replicateSupport(indelAssemblyPos1, indelRead, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO); // ensures this assembly is kept
        indelAssemblyNeg1.support().clear();
        indelAssemblyNeg1.addJunctionRead(indelRead);

        indelAssemblyNeg2.support().clear();
        replicateSupport(indelAssemblyNeg2, indelRead2Neg, ASSEMBLY_DEDUP_HIGH_SUPPORT_RATIO + 1); // would be kept but will be dropped

        dedupIndels.clear();

        existingAssemblies = Lists.newArrayList(indelAssemblyPos1);

        newAssemblies = Lists.newArrayList(indelAssemblyPos2);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(1, existingAssemblies.size());
        assertEquals(0, newAssemblies.size());
        assertTrue(dedupIndels.contains(indelAssemblyPos2));

        // ensure the neg assembly 2 is now dedup despite it having more support
        newAssemblies = Lists.newArrayList(indelAssemblyNeg1);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(1, existingAssemblies.size());
        assertEquals(1, newAssemblies.size());
        assertTrue(dedupIndels.contains(indelAssemblyPos2));

        existingAssemblies.add(indelAssemblyNeg1);

        newAssemblies = Lists.newArrayList(indelAssemblyNeg2);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(2, existingAssemblies.size());
        assertEquals(0, newAssemblies.size());
        assertTrue(dedupIndels.isEmpty());

        // test 4 - ensure that existing indels are only removed if the one being kept has a matching assembly
        existingAssemblies = Lists.newArrayList(indelAssemblyPos1, indelAssemblyNeg1);

        newAssemblies = Lists.newArrayList(indelAssemblyNeg2);

        dedupProximateAssemblies(existingAssemblies, newAssemblies, dedupIndels);

        assertEquals(1, existingAssemblies.size());
        assertTrue(existingAssemblies.contains(indelAssemblyPos1));
        assertEquals(1, newAssemblies.size());
    }

    private static void replicateSupport(final JunctionAssembly assembly, final Read read, int count)
    {
        for(int i = 0; i < count; ++i)
        {
            assembly.addJunctionRead(read);
        }
    }
}
