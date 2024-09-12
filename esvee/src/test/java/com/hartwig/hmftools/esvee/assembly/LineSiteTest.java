package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.TestUtils.makeCigarString;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.findLineExtensionEndIndex;
import static com.hartwig.hmftools.esvee.common.CommonUtils.withLineProximity;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_GAP;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_INDEL_MAX_OVERLAP;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.junit.Test;

public class LineSiteTest
{
    @Test
    public void testLineInsertionProximity()
    {
        assertFalse(withLineProximity(100, 100, FORWARD, FORWARD));

        assertTrue(withLineProximity(100, 100, FORWARD, REVERSE));
        assertTrue(withLineProximity(100, 100, REVERSE, FORWARD));

        // DEL within range
        int lowerPosition = 100;
        int upperPosition = lowerPosition + LINE_INDEL_MAX_GAP - 1;
        assertTrue(withLineProximity(lowerPosition, upperPosition, FORWARD, REVERSE));
        assertTrue(withLineProximity(upperPosition, lowerPosition, REVERSE, FORWARD));

        // too far for a DEL
        assertTrue(withLineProximity(lowerPosition, upperPosition + 1, FORWARD, REVERSE));
        assertTrue(withLineProximity(upperPosition + 1, lowerPosition, REVERSE, FORWARD));

        // DUP within range
        upperPosition = lowerPosition + LINE_INDEL_MAX_OVERLAP - 1;
        assertTrue(withLineProximity(lowerPosition, upperPosition, REVERSE, FORWARD));
        assertTrue(withLineProximity(upperPosition, lowerPosition, FORWARD, REVERSE));

        // too far for a DUP
        assertTrue(withLineProximity(lowerPosition, upperPosition + 1, REVERSE, FORWARD));
        assertTrue(withLineProximity(upperPosition + 1, lowerPosition, FORWARD, REVERSE));
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
    public void testLineAssemblyLinks()
    {
        String firstRefBases = REF_BASES_200.substring(0, 100);
        String secondRefBases = REF_BASES_200.substring(110, 200);

        String polyA = "AAAAAAAAAAAAAAAA";
        String polyT = Nucleotides.reverseComplementBases(polyA);

        // first a DEL
        Junction firstJunction = new Junction(CHR_1, 100, FORWARD);
        Junction secondJunction = new Junction(CHR_1, 110, REVERSE);

        String extraBases = "GTAGTGCTGTCGA";
        String firstExtBases = extraBases + polyA;
        String firstAssemblyBases = firstRefBases + firstExtBases;
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(firstAssemblyBases.length());

        JunctionAssembly firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, firstRefBases.length() - 1);

        String secondExtBases = polyA;
        String secondAssemblyBases = secondExtBases + secondRefBases;

        JunctionAssembly secondAssembly = new JunctionAssembly(secondJunction, secondAssemblyBases.getBytes(), baseQuals, secondExtBases.length());
        secondAssembly.markLineSequence();

        // order passed in doesn't matter
        AssemblyLink link = LineUtils.tryLineSequenceLink(firstAssembly, secondAssembly, false, false);
        assertNotNull(link);
        assertEquals(firstExtBases, link.insertedBases());
        assertEquals(firstAssembly, link.first());
        assertEquals(secondAssembly, link.second());
        assertEquals(DEL, link.svType());

        // same orientation
        firstExtBases = polyT + extraBases;
        firstAssemblyBases = firstExtBases + firstRefBases;

        firstJunction = new Junction(CHR_1, 100, REVERSE);
        firstAssembly = new JunctionAssembly(firstJunction, firstAssemblyBases.getBytes(), baseQuals, firstExtBases.length());

        link = LineUtils.tryLineSequenceLink(firstAssembly, secondAssembly, true, false);
        assertNotNull(link);
        assertEquals(firstExtBases, link.insertedBases());
        assertEquals(INV, link.svType());
    }
}