package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineExtensionEndIndex;
import static com.hartwig.hmftools.esvee.common.CommonUtils.withinLineProximity;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_GAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_OVERLAP;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.junit.Test;

public class LineSiteTest
{
    @Test
    public void testLineInsertionProximity()
    {
        assertFalse(withinLineProximity(100, 100, FORWARD, FORWARD));

        assertTrue(withinLineProximity(100, 100, FORWARD, REVERSE));
        assertTrue(withinLineProximity(100, 100, REVERSE, FORWARD));

        // DEL within range
        int lowerPosition = 100;
        int upperPosition = lowerPosition + LINE_INDEL_MAX_GAP - 1;
        assertTrue(withinLineProximity(lowerPosition, upperPosition, FORWARD, REVERSE));
        assertTrue(withinLineProximity(upperPosition, lowerPosition, REVERSE, FORWARD));

        // too far for a DEL
        assertTrue(withinLineProximity(lowerPosition, upperPosition + 1, FORWARD, REVERSE));
        assertTrue(withinLineProximity(upperPosition + 1, lowerPosition, REVERSE, FORWARD));

        // DUP within range
        upperPosition = lowerPosition + LINE_INDEL_MAX_OVERLAP - 1;
        assertTrue(withinLineProximity(lowerPosition, upperPosition, REVERSE, FORWARD));
        assertTrue(withinLineProximity(upperPosition, lowerPosition, FORWARD, REVERSE));

        // too far for a DUP
        assertTrue(withinLineProximity(lowerPosition, upperPosition + 1, REVERSE, FORWARD));
        assertTrue(withinLineProximity(upperPosition + 1, lowerPosition, FORWARD, REVERSE));
    }

    @Test
    public void testReadLineExtensionIndex()
    {
        String refBases = REF_BASES_200.substring(0, 50);
        String lineSequence = "AAAAAAAAAAAAAAAA";
        String extBases = lineSequence;
        String readBases = extBases + refBases;
        Read read = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, extBases.length(), 0));

        int lineIndex = findLineExtensionEndIndex(read, LINE_BASE_A, extBases.length() - 1, false);
        assertEquals(0, lineIndex);

        // test terminating after min requireed length
        extBases = "CCAAA" + lineSequence;
        readBases = extBases + refBases;
        read = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, extBases.length(), 0));

        lineIndex = findLineExtensionEndIndex(read, LINE_BASE_A, extBases.length() - 1, false);
        assertEquals(2, lineIndex);

        // test for forward junctions
        lineSequence = Nucleotides.reverseComplementBases("AAAAAAAAAAAAAAAA");
        extBases = lineSequence;
        readBases = refBases + extBases;
        read = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, 0, extBases.length()));

        lineIndex = findLineExtensionEndIndex(read, LINE_BASE_T, refBases.length(), true);
        assertEquals(readBases.length() - 1, lineIndex);

        extBases = lineSequence + "TTTAA";
        readBases = refBases + extBases;
        read = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, makeCigarString(readBases, 0, extBases.length()));

        lineIndex = findLineExtensionEndIndex(read, LINE_BASE_T, refBases.length(), true);
        assertEquals(readBases.length() - 3, lineIndex);
    }

    @Test
    public void testLineJunctionAssembly()
    {
        String polyAT = "TGTTTATTTTTTTTTTTT"; // allows for 2 non-bases

        String extraBases = "GTAGTGCTGT";

        Junction junction = new Junction(CHR_1, 100, FORWARD);

        // reads with poly-A tails of various lengths but all supporting the junction
        String refBases = REF_BASES_400.substring(51, 101);

        String extBases1 = polyAT + extraBases;
        String readBases1 = refBases + extBases1;
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases1, makeCigarString(readBases1, 0, extBases1.length()));
        Read read1b = cloneRead(read1, READ_ID_GENERATOR.nextId());

        String extBases2 = polyAT + "TT" + extraBases;
        String readBases2 = refBases.substring(1) + extBases2;
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 52, readBases2, makeCigarString(readBases2, 0, extBases2.length()));

        String extBases3 = polyAT + "T" + extraBases;
        String readBases3 = refBases + extBases3;
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases3, makeCigarString(readBases3, 0, extBases3.length()));

        // too short to be marked as line but supports the junction
        String extBases4 = polyAT.substring(0, 12);
        String readBases4 = refBases + extBases4;
        Read read4 = createRead(READ_ID_GENERATOR.nextId(), 51, readBases4, makeCigarString(readBases4, 0, extBases4.length()));

        List<Read> reads = List.of(read1, read1b, read2, read3, read4);
        reads.forEach(x -> ReadAdjustments.markLineSoftClips(x));

        assertFalse(read4.hasLineTail());

        JunctionAssembler junctionAssembler = new JunctionAssembler(junction);

        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());

        JunctionAssembly assembly = assemblies.get(0);

        assertTrue(assembly.hasLineSequence());
        assertEquals(5, assembly.supportCount());
        assertEquals(extBases1, assembly.formJunctionSequence()); // based on taking the median T length
    }

    @Test
    public void testLineAssemblyLinks()
    {
        String firstRefBases = REF_BASES_200.substring(0, 100);
        String secondRefBases = REF_BASES_200.substring(110, 200);

        String polyA = "AAAAAAAAAAAAAAAA";
        String polyT = Nucleotides.reverseComplementBases(polyA);

        // first a DEL with poly A at -ve junction
        Junction firstJunction = new Junction(CHR_1, 100, FORWARD);
        Junction secondJunction = new Junction(CHR_1, 110, REVERSE);

        String extraBases = "GTAGTGCTGTCGC";

        String firstExtBases = extraBases + polyA;
        String firstAssemblyBases = firstRefBases + firstExtBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(
                firstJunction, firstAssemblyBases.getBytes(), baseQuals, firstRefBases.length() - 1);

        String secondExtBases = polyA;
        String secondAssemblyBases = secondExtBases + secondRefBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(
                secondJunction, secondAssemblyBases.getBytes(), baseQuals, secondExtBases.length());
        secondAssembly.markLineSequence();

        AssemblyLink link = LineUtils.tryLineSequenceLink(firstAssembly, secondAssembly, false, false);
        assertNotNull(link);
        assertEquals(firstExtBases, link.insertedBases());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DEL, link.svType());

        // then a DEL with poly T at +ve junction
        firstExtBases = polyT;
        firstAssemblyBases = firstRefBases + firstExtBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, firstRefBases.length() - 1);
        firstAssembly.markLineSequence();

        secondExtBases = polyT + extraBases;
        secondAssemblyBases = secondExtBases + secondRefBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, secondExtBases.length());

        link = LineUtils.tryLineSequenceLink(firstAssembly, secondAssembly, false, false);
        assertNotNull(link);
        assertEquals(secondExtBases, link.insertedBases());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DEL, link.svType());

        // same orientation
        firstExtBases = polyT + extraBases;
        firstAssemblyBases = firstExtBases + firstRefBases;

        firstJunction = new Junction(CHR_1, 100, REVERSE);
        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, firstExtBases.length());

        secondExtBases = polyA;
        secondAssemblyBases = secondExtBases + secondRefBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, secondExtBases.length());
        secondAssembly.markLineSequence();

        link = LineUtils.tryLineSequenceLink(firstAssembly, secondAssembly, true, false);
        assertNotNull(link);
        String firstExtBasesRev = Nucleotides.reverseComplementBases(firstExtBases);
        assertEquals(firstExtBasesRev, link.insertedBases());
        assertEquals(INV, link.svType());

        // first assembly has the poly-T sequence
        String extraBases2 = "GTCGT";
        firstExtBases = polyT + extraBases2 + polyT;
        firstAssemblyBases = firstRefBases + firstExtBases;

        firstJunction = new Junction(CHR_1, 120, FORWARD);
        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, firstRefBases.length() - 1);
        firstAssembly.markLineSequence();

        secondExtBases = polyT + extraBases2 + polyT + extraBases;
        secondAssemblyBases = secondExtBases + secondRefBases;

        secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, secondExtBases.length());

        link = LineUtils.tryLineSequenceLink(firstAssembly, secondAssembly, false, false);
        assertNotNull(link);
        assertEquals(firstExtBases + extraBases, link.insertedBases());
        assertEquals(DUP, link.svType());
    }
}