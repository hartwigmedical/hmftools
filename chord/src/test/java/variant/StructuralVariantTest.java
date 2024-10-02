package variant;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.chord.variant.StructuralVariant;

import org.junit.Test;

public class StructuralVariantTest
{
    @Test
    public void canParseDelAndDupLength()
    {
        // Check both start and end since bracket orientations are reversed for DEL and DUP ends (higher position in genome)

        // Deletion
        StructuralVariant delStart = new StructuralVariant("1", 207981231, "A", "A[1:208014820[");
        StructuralVariant delEnd   = new StructuralVariant("1", 208014820, "G", "]1:207981231]G");

        assertEquals(DEL, delStart.Type);
        assertEquals(DEL, delEnd.Type);
        assertEquals(33589, delStart.Length.intValue());
        assertEquals(33589, delStart.Length.intValue());

        // Duplication
        StructuralVariant dupStart = new StructuralVariant("4", 66212039, "T", "T[4:66211957[");
        StructuralVariant dupEnd   = new StructuralVariant("4", 66211957, "G", "]4:66212039]G");

        assertEquals(DUP, dupStart.Type);
        assertEquals(DUP, dupEnd.Type);
        assertEquals(82, dupStart.Length.intValue());
        assertEquals(82, dupEnd.Length.intValue());
    }

    @Test
    public void canParseInvAndBnd()
    {
        // Inversion
        StructuralVariant invStart = new StructuralVariant("3", 24565108, "A", "[3:24566180[A");
        assertEquals(INV, invStart.Type);
        assertEquals(1072, invStart.Length.intValue());

        // Translocation
        StructuralVariant bndStart = new StructuralVariant("1", 87337011, "T", "T]10:36119127]");
        assertEquals(BND, bndStart.Type);
        assertNull(bndStart.Length);
    }

    @Test
    public void canParseUnusedSvTypes()
    {
        // Insertion
        StructuralVariant insStart = new StructuralVariant("3", 60204585, "C", "CTACATACATACATACATACATACATACATACA[3:60204586[");
        assertEquals(INS, insStart.Type);
        assertEquals(1, insStart.Length.intValue());

        // Single breakend
        StructuralVariant sglStart = new StructuralVariant("2", 845380, "C", "CAGGTCATGGGACATTTCCCTTGAGTGTTCACAGGTGCACAGAAACACGGGCCATGGGACATTTCCCACGA.");
        assertEquals(SGL, sglStart.Type);
        assertNull(sglStart.Length);
    }
}
