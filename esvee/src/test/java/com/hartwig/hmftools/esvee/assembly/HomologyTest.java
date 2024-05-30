package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_GENOME;
import static com.hartwig.hmftools.esvee.alignment.HomologyData.determineHomology;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAlignment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.HomologyData;

import org.junit.Test;

public class HomologyTest
{
    @Test
    public void testHomology()
    {
        AlignData alignmentStart = createAlignment(CHR_1, 100, 150, 0, 50, "51M");
        AlignData alignmentEnd = createAlignment(CHR_1, 100, 150, 51, 100, "50M");

        String assemblyOverlap = "";

        assertNull(determineHomology(assemblyOverlap, alignmentStart, alignmentEnd, REF_GENOME));

        String basesStart = "TTCTTCTTCTC";
        String basesEnd = basesStart;
        assemblyOverlap = basesStart;

        // test 1: exact match
        HomologyData homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertNotNull(homology);
        assertEquals(basesStart, homology.Homology);
        assertEquals(-6, homology.ExactStart);
        assertEquals(5, homology.ExactEnd);
        assertEquals(-6, homology.InexactStart);
        assertEquals(5, homology.InexactEnd);

        // test 2: first base now no longer matches
        basesStart = "GTCTTCTTCTC";
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("", homology.Homology);
        assertEquals(0, homology.ExactStart);
        assertEquals(0, homology.ExactEnd);
        assertEquals(0, homology.InexactStart);
        assertEquals(11, homology.InexactEnd);

        // test 3: first base matches, range of lowest mismatches is 0-1
        basesStart = "TGCATCTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = assemblyOverlap;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("T", homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(0, homology.ExactEnd);
        assertEquals(-1, homology.InexactStart);
        assertEquals(10, homology.InexactEnd);

        // test 4: second base has mismatch
        basesStart = "TTGATCTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = assemblyOverlap;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("TT", homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(1, homology.ExactEnd);
        assertEquals(-1, homology.InexactStart);
        assertEquals(10, homology.InexactEnd);

        // test 5: min mismatches in range 4-6
        basesStart = "TTCATCGTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = "TTCTTCTTCTC";
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("TTC", homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(1, homology.ExactEnd);
        assertEquals(-5, homology.InexactStart);
        assertEquals(6, homology.InexactEnd);

        // test 6: ref sequences match but have the same difference from the assembly
        basesStart = "TTCATGTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = basesStart;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals(basesStart, homology.Homology);
        assertEquals(-6, homology.ExactStart);
        assertEquals(5, homology.ExactEnd);
        assertEquals(-6, homology.InexactStart);
        assertEquals(5, homology.InexactEnd);

        // test 7:
        basesStart = "TTCATGTTCTC";
        assemblyOverlap = "TTCATCTTCTC";
        basesEnd = "TTCATATTCTC";
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals("TTCAT", homology.Homology);
        assertEquals(-6, homology.ExactStart);
        assertEquals(5, homology.ExactEnd);
        assertEquals(-6, homology.InexactStart);
        assertEquals(5, homology.InexactEnd);

        // test 8: with an even number of bases
        basesStart = "AA";
        assemblyOverlap = basesStart;
        basesEnd = basesStart;
        homology = determineHomology(assemblyOverlap, basesStart, basesEnd, basesStart.length());
        assertEquals(basesStart, homology.Homology);
        assertEquals(-1, homology.ExactStart);
        assertEquals(1, homology.ExactEnd);
        assertEquals(-1, homology.InexactStart);
        assertEquals(1, homology.InexactEnd);
    }

    @Test
    public void testPositionShift()
    {
        HomologyData homology = new HomologyData("", -2, 1, -2, 1);
        assertEquals(-1, homology.positionAdjustment(FORWARD, true));
        assertEquals(2, homology.positionAdjustment(REVERSE, false));

        // the same-orientation cases
        assertEquals(1, homology.positionAdjustment(REVERSE, true));
        assertEquals(-2, homology.positionAdjustment(FORWARD, false));
    }


}
