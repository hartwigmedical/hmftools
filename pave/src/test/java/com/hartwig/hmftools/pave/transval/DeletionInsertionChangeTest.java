package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class DeletionInsertionChangeTest extends TransvalTestBase
{
    @Test
    public void prefixLengthTest()
    {
        assertEquals(0, DeletionInsertionChange.lengthOfCommonPrefix("", ""));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonPrefix("", "W"));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonPrefix("", "Whatever"));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonPrefix("W", ""));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonPrefix("Whatever", ""));
        assertEquals(1, DeletionInsertionChange.lengthOfCommonPrefix("Whatever", "Water"));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonPrefix("Whatever", "water"));
        assertEquals(2, DeletionInsertionChange.lengthOfCommonPrefix("Whatever", "Which"));
        assertEquals(4, DeletionInsertionChange.lengthOfCommonPrefix("Whatever", "What"));
        assertEquals(4, DeletionInsertionChange.lengthOfCommonPrefix("What", "Whatever"));
        assertEquals(4, DeletionInsertionChange.lengthOfCommonPrefix("What", "What"));
    }

    @Test
    public void suffixLength()
    {
        assertEquals(0, DeletionInsertionChange.lengthOfCommonSuffix("", ""));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonSuffix("", "r"));
        assertEquals(0, DeletionInsertionChange.lengthOfCommonSuffix("r", ""));
        assertEquals(1, DeletionInsertionChange.lengthOfCommonSuffix("r", "er"));
        assertEquals(1, DeletionInsertionChange.lengthOfCommonSuffix("er", "r"));
        assertEquals(2, DeletionInsertionChange.lengthOfCommonSuffix("er", "er"));
        assertEquals(4, DeletionInsertionChange.lengthOfCommonSuffix("Whoever", "Whatever"));
    }

    @Test
    public void noMatches()
    {
        DeletionInsertionChange change = change("AAAAAA", "CCC");
        assertEquals("AAAAAA", change.deleted());
        assertEquals("CCC", change.inserted());
        assertEquals(0, change.positionOfDeletion());

        change = change("AAAAA", "CCC");
        assertEquals("AAAAA", change.deleted());
        assertEquals("CCC", change.inserted());
        assertEquals(0, change.positionOfDeletion());

        change = change("AAAA", "CCC");
        assertEquals("AAAA", change.deleted());
        assertEquals("CCC", change.inserted());
        assertEquals(0, change.positionOfDeletion());

        change = change("AACCCAA", "CCC");
        assertEquals("AACCCAA", change.deleted());
        assertEquals("CCC", change.inserted());
        assertEquals(0, change.positionOfDeletion());

        change = change("ACCCA", "CCC");
        assertEquals("ACCCA", change.deleted());
        assertEquals("CCC", change.inserted());
        assertEquals(0, change.positionOfDeletion());
    }

    @Test
    public void rightSidesMatch()
    {
        DeletionInsertionChange change = change("GCGGAGAACTGG", "ATG");
        assertEquals("GCGGAGAACTG", change.deleted());
        assertEquals("AT", change.inserted());
        assertEquals(0, change.positionOfDeletion());

        change = change("TTAAGAGAAGCA", "CCA");
        assertEquals("TTAAGAGAAG", change.deleted());
        assertEquals("C", change.inserted());
        assertEquals(0, change.positionOfDeletion());
    }

    @Test
    public void leftSidesMatch()
    {
        // VHL A5_W8 => DT
        DeletionInsertionChange change = change("GCGGAGAACTGG", "GATACC");
        assertEquals("CGGAGAACTGG", change.deleted());
        assertEquals("ATACC", change.inserted());
        assertEquals(1, change.positionOfDeletion());

        // VHL A5_W8 => AT
        change = change("GCGGAGAACTGG", "GCTACC");
        assertEquals("GGAGAACTGG", change.deleted());
        assertEquals("TACC", change.inserted());
        assertEquals(2, change.positionOfDeletion());
    }

    @Test
    public void leftMatchIsLongest()
    {
        // VHL A5_W8 => AT
        DeletionInsertionChange change = change("GCGGAGAACTGG", "GCTACG");
        assertEquals("GGAGAACTG", change.deleted());
        assertEquals("TAC", change.inserted());
        assertEquals(2, change.positionOfDeletion());
    }

    @Test
    public void chooseLeftMost()
    {
        DeletionInsertionChange change = change("GAATAT", "GAT");
        assertEquals("AAT", change.deleted());
        assertEquals("", change.inserted());
        assertEquals(1, change.positionOfDeletion());
    }

    @Test
    public void rightMatchIsLongest()
    {
        // VHL A5_W8 => DR
        DeletionInsertionChange change = change("GCGGAGAACTGG", "GATAGG");
        assertEquals("CGGAGAACT", change.deleted());
        assertEquals("ATA", change.inserted());
        assertEquals(1, change.positionOfDeletion());
    }

    @Test
    public void leftAndRightMatchesAreEqual()
    {
        // VHL A5_W8 => DT
        DeletionInsertionChange change = change("GCGGAGAACTGG", "GATACG");
        assertEquals("CGGAGAACTG", change.deleted());
        assertEquals("ATAC", change.inserted());
        assertEquals(1, change.positionOfDeletion());

        // VHL A5_W8 => AW
        change = change("GCGGAGAACTGG", "GCTAGG");
        assertEquals("GGAGAACT", change.deleted());
        assertEquals("TA", change.inserted());
        assertEquals(2, change.positionOfDeletion());
    }

    private DeletionInsertionChange change(String ref, String alt)
    {
        return new DeletionInsertionChange(ref, alt);
    }
}
