package com.hartwig.hmftools.pavereverse.variants;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

public class CodonChangeTest
{
    @Test(expected = IllegalArgumentException.class)
    public void cannotCompareVariantsWithDifferingReferenceCodons()
    {
        cv("AAA", "TTT").compareTo(cv("AAG", "TTT"));
    }

    @Test
    public void compareTest()
    {
        assertTrue(cv("AAA", "AAA").compareTo(cv("AAA", "AAC")) < 0);
        assertTrue(cv("AAT", "AAT").compareTo(cv("AAT", "AAC")) < 0);
        assertEquals(0, cv("AAT", "AAC").compareTo(cv("AAT", "AAC")));

        assertTrue(cv("AAT", "AAC").compareTo(cv("AAT", "AAA")) > 0);
        assertTrue(cv("AAT", "AAC").compareTo(cv("AAT", "AAG")) < 0);

        assertTrue(cv("AAT", "AAC").compareTo(cv("AAT", "ACG")) < 0);
        assertTrue(cv("TAT", "AAC").compareTo(cv("TAT", "AAG")) < 0);
    }

    @Test
    public void editDistanceTest()
    {
        assertEquals(0, cv("AAA", "AAA").editDistance());
        assertEquals(1, cv("AAA", "AAG").editDistance());
        assertEquals(1, cv("AAA", "AAT").editDistance());
        assertEquals(1, cv("AAA", "AGA").editDistance());
        assertEquals(1, cv("AAA", "CAA").editDistance());
        assertEquals(2, cv("AAA", "CAG").editDistance());
        assertEquals(2, cv("AAA", "ATG").editDistance());
        assertEquals(2, cv("AAA", "TTA").editDistance());
        assertEquals(3, cv("AAA", "TTG").editDistance());
        assertEquals(3, cv("AAA", "TTC").editDistance());
    }

    @Test
    public void differenceStringsTest()
    {
        assertEquals(p("", ""), cv("AAA", "AAA").differenceStrings());

        assertEquals(p("A", "T"), cv("AAA", "AAT").differenceStrings());
        assertEquals(p("A", "T"), cv("AAA", "ATA").differenceStrings());
        assertEquals(p("A", "T"), cv("AAA", "TAA").differenceStrings());

        assertEquals(p("GC", "TT"), cv("AGC", "ATT").differenceStrings());
        assertEquals(p("GC", "TT"), cv("GCA", "TTA").differenceStrings());
        assertEquals(p("GAC", "TAT"), cv("GAC", "TAT").differenceStrings());

        assertEquals(p("GCC", "TTA"), cv("GCC", "TTA").differenceStrings());
    }

    @Test
    public void testEquals()
    {
        assertEquals(cv("AAA", "AAG"), cv("AAA", "AAG"));
        assertNotEquals(cv("AAA", "AAC"), cv("AAA", "AAG"));
        assertNotEquals(cv("AAC", "AAG"), cv("AAA", "AAG"));
    }
    
    @Test
    public void testHashCode()
    {
        assertEquals(cv("AAA", "AAG").hashCode(), cv("AAA", "AAG").hashCode());
    }
    
    private CodonChange cv(String s, String t)
    {
        return new CodonChange(s, t);
    }

    private Pair<String,String> p(String s, String t)
    {
        return Pair.of(s, t);
    }
}
