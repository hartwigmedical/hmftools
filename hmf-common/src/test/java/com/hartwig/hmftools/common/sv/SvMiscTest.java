package com.hartwig.hmftools.common.sv;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_NEG_CHAR;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_POS_CHAR;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formSingleAltString;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.region.Orientation;

import org.junit.Test;

public class SvMiscTest
{
    @Test
    public void testOrientation()
    {
        assertEquals(FORWARD, REVERSE.opposite());
        assertEquals(REVERSE, FORWARD.opposite());

        assertTrue(FORWARD.equalsByte(ORIENT_FWD));
        assertTrue(REVERSE.equalsByte(ORIENT_REV));

        assertEquals(FORWARD, Orientation.fromByte(ORIENT_FWD));
        assertEquals(REVERSE, Orientation.fromByte(ORIENT_REV));

        assertEquals(FORWARD, Orientation.fromChar(ORIENT_POS_CHAR));
        assertEquals(REVERSE, Orientation.fromChar(ORIENT_NEG_CHAR));

        assertEquals(FORWARD, Orientation.fromByteStr(String.valueOf(ORIENT_FWD)));
        assertEquals(REVERSE, Orientation.fromByteStr(String.valueOf(ORIENT_REV)));
    }

    @Test
    public void testSvBreakendVcfFormation()
    {
        // pos orientation for DEL and DUP: AGAGATTATACTTTGTGTA[10:89712341[
        // pos orientation for INV: G]3:26664499]

        // neg orientation for DEL and DUP: ]10:89700299]GAGATTATACTTTGTGTAA
        // neg orientation for INV: [3:24566181[C

        String alt = "A";
        String insert = "GT";
        String otherChromosome = "2";
        int otherPosition = 10000;
        Orientation orientation = FORWARD;
        Orientation otherOrientation = REVERSE;

        String altStr = formPairedAltString(alt, insert, otherChromosome, otherPosition, orientation, otherOrientation);

        assertEquals("AGT[2:10000[", altStr);

        VariantAltInsertCoords altCoords = fromRefAlt(altStr, alt);

        assertEquals(alt, altCoords.Alt);
        assertEquals(insert, altCoords.InsertSequence);
        assertEquals(orientation, altCoords.Orient);
        assertEquals(otherChromosome, altCoords.OtherChromsome);
        assertEquals(otherPosition, altCoords.OtherPosition);
        assertEquals(otherOrientation, altCoords.OtherOrient);

        orientation = REVERSE;
        otherOrientation = FORWARD;

        altStr = formPairedAltString(alt, insert, otherChromosome, otherPosition, orientation, otherOrientation);

        assertEquals("]2:10000]GTA", altStr);

        altCoords = fromRefAlt(altStr, alt);

        assertEquals(alt, altCoords.Alt);
        assertEquals(insert, altCoords.InsertSequence);
        assertEquals(orientation, altCoords.Orient);
        assertEquals(otherChromosome, altCoords.OtherChromsome);
        assertEquals(otherPosition, altCoords.OtherPosition);
        assertEquals(otherOrientation, altCoords.OtherOrient);

        // same-orientation breakends
        orientation = FORWARD;
        otherOrientation = FORWARD;

        altStr = formPairedAltString(alt, insert, otherChromosome, otherPosition, orientation, otherOrientation);

        assertEquals("AGT]2:10000]", altStr);

        altCoords = fromRefAlt(altStr, alt);

        assertEquals(alt, altCoords.Alt);
        assertEquals(insert, altCoords.InsertSequence);
        assertEquals(orientation, altCoords.Orient);
        assertEquals(otherChromosome, altCoords.OtherChromsome);
        assertEquals(otherPosition, altCoords.OtherPosition);
        assertEquals(otherOrientation, altCoords.OtherOrient);

        orientation = REVERSE;
        otherOrientation = REVERSE;

        altStr = formPairedAltString(alt, insert, otherChromosome, otherPosition, orientation, otherOrientation);

        assertEquals("[2:10000[GTA", altStr);

        altCoords = fromRefAlt(altStr, alt);

        assertEquals(alt, altCoords.Alt);
        assertEquals(insert, altCoords.InsertSequence);
        assertEquals(orientation, altCoords.Orient);
        assertEquals(otherChromosome, altCoords.OtherChromsome);
        assertEquals(otherPosition, altCoords.OtherPosition);
        assertEquals(otherOrientation, altCoords.OtherOrient);

        // single breakends
        orientation = FORWARD;

        altStr = formSingleAltString(alt, insert, orientation);

        assertEquals("AGT.", altStr);

        altCoords = fromRefAlt(altStr, alt);

        assertEquals(alt, altCoords.Alt);
        assertEquals(insert, altCoords.InsertSequence);
        assertEquals(orientation, altCoords.Orient);

        orientation = REVERSE;

        altStr = formSingleAltString(alt, insert, orientation);

        assertEquals(".GTA", altStr);

        altCoords = fromRefAlt(altStr, alt);

        assertEquals(alt, altCoords.Alt);
        assertEquals(insert, altCoords.InsertSequence);
        assertEquals(orientation, altCoords.Orient);
    }

}
