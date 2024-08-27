package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_A;
import static com.hartwig.hmftools.common.sv.LineElements.LINE_BASE_T;
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
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.esvee.assembly.read.Read;

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
}
