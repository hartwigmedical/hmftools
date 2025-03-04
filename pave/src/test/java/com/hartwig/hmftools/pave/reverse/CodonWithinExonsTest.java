package com.hartwig.hmftools.pave.reverse;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import org.junit.Test;

public class CodonWithinExonsTest extends ReversePaveTestBase
{
    @Test
    public void possibleVariantsGivingTest()
    {
        CodonWithinExons codon = new CodonWithinExons(bs(100, "ATC", true));
        Set<CodonChange> changes = codon.possibleVariantsGiving(aa("M"));
        assertEquals(1, changes.size());
        CodonChange change = changes.iterator().next();
        assertEquals("ATC",change.ReferenceCodon);
        assertEquals("ATG",change.AlternateCodon);
        assertEquals(1,change.editDistance());

        changes = codon.possibleVariantsGiving(aa("Y"));
        assertEquals(2, changes.size());
        assertTrue(changes.contains(new CodonChange("ATC", "TAC")));
        assertTrue(changes.contains(new CodonChange("ATC", "TAT")));

        codon = new CodonWithSuffixInNextExon(bs(100, "AT", true), "C");
        changes = codon.possibleVariantsGiving(aa("Y"));
        assertEquals(1, changes.size());
        assertTrue(changes.contains(new CodonChange("ATC", "TAC")));

        codon = new CodonWithPrefixInPreviousExon("T", bs(100, "AT", true));
        changes = codon.possibleVariantsGiving(aa("L"));
        assertEquals(2, changes.size());
        assertTrue(changes.contains(new CodonChange("TAT", "TTA")));
        assertTrue(changes.contains(new CodonChange("TAT", "TTG")));

        changes = codon.possibleVariantsGiving(aa("E"));
        assertTrue(changes.isEmpty());
    }

    @Test
    public void possibleVariantsGivingStopTest()
    {
        CodonWithinExons codon = new CodonWithinExons(bs(100, "ATC", true));
        Set<CodonChange> changes = codon.possibleVariantsGivingStop();
        assertEquals(3, changes.size());
        assertTrue(changes.contains(new CodonChange("ATC", "TAA")));
        assertTrue(changes.contains(new CodonChange("ATC", "TAG")));
        assertTrue(changes.contains(new CodonChange("ATC", "TGA")));

        codon = new CodonWithPrefixInPreviousExon("G", bs(100, "AT", true));
        changes = codon.possibleVariantsGivingStop();
        assertTrue(changes.isEmpty());
    }
}
