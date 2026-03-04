package com.hartwig.hmftools.finding;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class HlaAlleleFactoryTest
{
    @Test
    public void testRegEx()
    {
        var matcher = HlaAlleleFactory.matchHlaRegEx("A*01:01");

        //throw IllegalStateException("Can't extract HLA gene, alleleGroup and hlaProtein from ${allele.allele()}")
        assertEquals("A", matcher.group("gene"));
        assertEquals("01", matcher.group("alleleGroup"));
        assertEquals("01", matcher.group("hlaProtein"));
    }
}