package com.hartwig.hmftools.finding;

import static com.hartwig.hmftools.finding.HlaAlleleFactory.HLA_REGEX;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class HlaAlleleFactoryTest
{
    @Test
    public void testRegEx()
    {
        var matcher = HLA_REGEX.matcher("A*01:01");
        assertTrue(matcher.matches());

        //throw IllegalStateException("Can't extract HLA gene, alleleGroup and hlaProtein from ${allele.allele()}")
        assertEquals("A", matcher.group("gene"));
        assertEquals("01", matcher.group("alleleGroup"));
        assertEquals("01", matcher.group("hlaProtein"));
    }
}