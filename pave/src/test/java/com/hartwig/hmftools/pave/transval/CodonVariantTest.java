package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class CodonVariantTest
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
    
    private CodonVariant cv(String s, String t)
    {
        return new CodonVariant(s, t);
    }
}
