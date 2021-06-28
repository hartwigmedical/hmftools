package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.hla.HlaAllele.dedup;
import static com.hartwig.hmftools.lilac.seq.HlaSequenceFile.asSixDigit;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class HlaAlleleTest
{
    @Test
    public void testDecode()
    {
        assertContig("A*01:01:01:01", "A*01:01:01:01");
        assertContig("A*01:01:01", "A*01:01:01");
        assertContig("A*26:22", "A*26:22");
        assertContig("A*26", "A*26");
    }

    @Test
    public void testReduce()
    {
        HlaAllele eightDigit = HlaAllele.fromString("A*01:01:01:01");
        assertMatch(HlaAllele.fromString("A*01:01:01"), asSixDigit(eightDigit));
        assertMatch(HlaAllele.fromString("A*01:01"), eightDigit.asFourDigit());
        assertMatch(HlaAllele.fromString("A*01"), eightDigit.asAlleleGroup());
    }

    @Test
    public void testDedupAlleles()
    {
        List<HlaAllele> alleles = Lists.newArrayList();
        HlaAllele allele1 = HlaAllele.fromString("A*01:01:01");
        HlaAllele allele2 = HlaAllele.fromString("A*01:02:01");
        HlaAllele allele3 = HlaAllele.fromString("A*01:03:01");
        alleles.add(allele1);
        alleles.add(allele1);
        alleles.add(allele2);
        alleles.add(allele2);
        alleles.add(allele3);
        List<HlaAllele> deduped = dedup(alleles);
        assertEquals(3, deduped.size());
        assertTrue(deduped.contains(allele1));
        assertTrue(deduped.contains(allele2));
        assertTrue(deduped.contains(allele3));
    }

    @Test
    public void testAlleleCache()
    {
        HlaAlleleCache cache = new HlaAlleleCache();
        HlaAllele allele1 = cache.request("A*01:01:01");
        HlaAllele allele2 = cache.request("A*01:01:02");
        HlaAllele allele3 = cache.request("A*01:02:01");
        HlaAllele allele4 = cache.request("A*01:02:02");
        assertEquals(4, cache.alleleCount());
        assertEquals(2, cache.fourDigitCount());
        assertEquals(1, cache.groupCount());

        HlaAllele fourDigit1 = cache.requestFourDigit("A*01:01");
        assertEquals(fourDigit1, allele1.asFourDigit());
        assertEquals(fourDigit1, allele2.asFourDigit());
        assertEquals(allele3.asFourDigit(), allele4.asFourDigit());
        assertEquals(allele1.asAlleleGroup(), allele4.asAlleleGroup());

        HlaAllele group1 = cache.requestGroup(HlaAllele.fromString("A*01"));
        HlaAllele fourDigit2 = cache.requestFourDigit("A*01:02");
        assertEquals(fourDigit2, allele3.asFourDigit());

        assertEquals(group1, allele1.asAlleleGroup());
        assertEquals(group1, fourDigit1.asAlleleGroup());
        assertEquals(2, cache.fourDigitCount());
        assertEquals(1, cache.groupCount());
    }

    @Test
    public void testAlleleMaps()
    {
        Map<HlaAllele,Integer> map = Maps.newHashMap();

        HlaAllele allele1 = HlaAllele.fromString("A*01:01");
        HlaAllele allele2 = HlaAllele.fromString("A*01:01");
        HlaAllele allele3 = HlaAllele.fromString("A*01:01:01");
        HlaAllele allele4 = HlaAllele.fromString("A*01:02");

        map.put(allele1, 1);
        map.put(allele2, 2);
        map.put(allele3, 3);
        map.put(allele4, 4);

        assertEquals(map.size(), 3);

        Integer val1 = map.get(allele1);
        Integer val2 = map.get(allele2);

        assertEquals(val1, val2);

        map.remove(allele1);
        assertEquals(map.size(), 2);

        map.remove(allele3);
        assertEquals(map.size(), 1);

        map.remove(allele4);
        assertTrue(map.isEmpty());
    }

    private void assertContig(String expected, String contig)
    {
        assertEquals(expected, HlaAllele.fromString(contig).toString());
    }

    private void assertMatch(HlaAllele expected, HlaAllele contig)
    {
        assertTrue(expected.toString().equals(contig.toString()));
    }

}
