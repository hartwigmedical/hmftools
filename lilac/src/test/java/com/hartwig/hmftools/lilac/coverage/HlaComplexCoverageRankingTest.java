package com.hartwig.hmftools.lilac.coverage;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class HlaComplexCoverageRankingTest
{
    @Test
    public void testHomozygousCounts()
    {
        HlaComplexCoverage coverage = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0)));

        assertEquals(3, coverage.homozygousCount());

        coverage = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:03"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0)));

        assertEquals(1, coverage.homozygousCount());

        coverage = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:03"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*03:02"), 10, 0, 0)));

        assertEquals(0, coverage.homozygousCount());
    }

    @Test
    public void testCoverageRankSorting()
    {
        HlaComplexCoverage coverage1 = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:201"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*01:01"), 10, 0, 0)));

        coverage1.setScore(34);

        // middle 2 are sorted by their alleles numerically
        HlaComplexCoverage coverage2 = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*02:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*02:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*02:01"), 10, 0, 0)));

        coverage2.setScore(33);

        HlaComplexCoverage coverage3 = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*03:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*03:01"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0)));

        coverage3.setScore(35);

        HlaComplexCoverage coverage4 = HlaComplexCoverage.create(Lists.newArrayList(
                new HlaAlleleCoverage(HlaAllele.fromString("A*01:78"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0, 0),
                new HlaAlleleCoverage(HlaAllele.fromString("C*01:02"), 10, 0, 0)));

        coverage4.setScore(34);

        List<HlaComplexCoverage> complexCoverages = Lists.newArrayList(coverage1, coverage2, coverage3, coverage4);

        Collections.sort(complexCoverages, new HlaComplexCoverageRanking.ComplexCoverageSorter());
        assertEquals(coverage3, complexCoverages.get(0));
        assertEquals(coverage4, complexCoverages.get(1));
        assertEquals(coverage1, complexCoverages.get(2));
        assertEquals(coverage2, complexCoverages.get(3));
    }


    /*

    private HlaComplexCoverageRankingOld victim = new HlaComplexCoverageRankingOld(
            3, Lists.newArrayList(), Lists.newArrayList(), Lists.newArrayList());

    private HlaAlleleCoverage a1 = new HlaAlleleCoverage(HlaAllele.fromString("A*01:01"), 10, 0.0, 0.0);
    private HlaAlleleCoverage a2 = new HlaAlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0.0, 0.0);
    private HlaAlleleCoverage a3 = new HlaAlleleCoverage(HlaAllele.fromString("A*01:03"), 10, 0.0, 0.0);
    private HlaAlleleCoverage b1 = new HlaAlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0.0, 0.0);
    private HlaAlleleCoverage b2 = new HlaAlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0.0, 0.0);
    private HlaAlleleCoverage b3 = new HlaAlleleCoverage(HlaAllele.fromString("B*01:03"), 10, 0.0, 0.0);
    private HlaAlleleCoverage c1 = new HlaAlleleCoverage(HlaAllele.fromString("C*01:01"), 10, 0.0, 0.0);
    private HlaAlleleCoverage c2 = new HlaAlleleCoverage(HlaAllele.fromString("C*01:02"), 10, 0.0, 0.0);
    private HlaAlleleCoverage c3 = new HlaAlleleCoverage(HlaAllele.fromString("C*01:03"), 10, 0.0, 0.0);

    @Test
    public void testMostCommon()
    {
        HlaComplexCoverage common = HlaComplexCoverage.create(Lists.newArrayList(a1, a3, b1, b3, c1, c2));
        List<HlaComplexCoverage> lessCommon = Lists.newArrayList();

        lessCommon.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c2)));
        lessCommon.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a3, b1, b2, c1, c2)));
        lessCommon.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b3, c1, c2)));
        lessCommon.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c3)));

        for(HlaComplexCoverage lower : lessCommon)
        {
            List<HlaComplexCoverage> complexes = Lists.newArrayList(common, lower); //.shuffled()

            HlaComplexCoverage winner = new HlaComplexCoverageRankingOld(
                    3, Lists.newArrayList(a3.Allele, b3.Allele, c3.Allele),
                    Lists.newArrayList(), Lists.newArrayList()).candidateRanking(complexes).get(0);
            assertEquals(common, winner);
        }
    }

    @Test
    public void testLeastRecovered()
    {
        List<HlaAllele> common = Lists.newArrayList(a1.Allele, a2.Allele, a3.Allele);
        List<HlaAllele> recovered = Lists.newArrayList(a1.Allele);

        HlaComplexCoverage victim1 = HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c2));
        HlaComplexCoverage victim2 = HlaComplexCoverage.create(Lists.newArrayList(a1, a3, b1, b2, c1, c2));
        HlaComplexCoverage victim3 = HlaComplexCoverage.create(Lists.newArrayList(a2, a3, b1, b2, c1, c2));
        List<HlaComplexCoverage> complexes = Lists.newArrayList(victim1, victim2, victim3); // .shuffled()

        List<HlaComplexCoverage> ranked = new HlaComplexCoverageRankingOld(
                3, common, recovered, Lists.newArrayList()).candidateRanking(complexes);
        HlaComplexCoverage winner = ranked.get(0);
        assertEquals(victim3, winner);
    }

    @Test
    public void testHighestAllele()
    {
        HlaComplexCoverage highest = HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c2));
        List<HlaComplexCoverage> lowerList = Lists.newArrayList();

        lowerList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a3, b1, b2, c1, c2)));
        lowerList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b3, c1, c2)));
        lowerList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c3)));

        for (HlaComplexCoverage lower : lowerList)
        {
            List<HlaComplexCoverage> complexes = Lists.newArrayList(highest, lower); //.shuffled()
            HlaComplexCoverage winner = victim.candidateRanking(complexes).get(0);
            assertEquals(highest, winner);
        }
    }

    @Test
    public void testHomozygousWinsIfAllElseEqual()
    {
        HlaComplexCoverage het = HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c2));
        List<HlaComplexCoverage> homList = Lists.newArrayList();

        homList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a1, b1, b2, c1, c2)));
        homList.add(HlaComplexCoverage.create(Lists.newArrayList(a2, a2, b1, b2, c1, c2)));
        homList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b1, c1, c2)));
        homList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b2, b2, c1, c2)));
        homList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c1)));
        homList.add(HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c2, c2)));

        for (HlaComplexCoverage hom : homList)
        {
            List<HlaComplexCoverage> complexes = Lists.newArrayList(het, hom); // .shuffled()
            HlaComplexCoverage winner = victim.candidateRanking(complexes).get(0);
            assertEquals(hom, winner);
        }

        HlaComplexCoverage homA = HlaComplexCoverage.create(Lists.newArrayList(a1, a1, b1, b2, c1, c2));
        HlaComplexCoverage homAAndB = HlaComplexCoverage.create(Lists.newArrayList(a1, a1, b2, b2, c1, c2));
        List<HlaComplexCoverage> complexes = Lists.newArrayList(homA, homAAndB); // .shuffled();
        assertEquals(homAAndB, victim.candidateRanking(complexes).get(0));
    }

    @Test
    public void testWildCoverage()
    {
        HlaComplexCoverage het = HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c2));
        HlaAlleleCoverage c2Copy = new HlaAlleleCoverage(c2.Allele, 9, c2.SharedCoverage, 1);
        HlaComplexCoverage hetWithOneWild = HlaComplexCoverage.create(Lists.newArrayList(a1, a2, b1, b2, c1, c2Copy));

        HlaComplexCoverage hom = HlaComplexCoverage.create(Lists.newArrayList(a1, a1, b1, b1, c1, c1));
        HlaAlleleCoverage c1Copy = new HlaAlleleCoverage(c1.Allele, 9, c1.SharedCoverage, 1);
        HlaComplexCoverage homWithOneWild = HlaComplexCoverage.create(Lists.newArrayList(a1, a1, b1, b1, c1, c1Copy));
        HlaAlleleCoverage b1Copy = new HlaAlleleCoverage(b1.Allele, 9, b1.SharedCoverage, 1);
        HlaComplexCoverage homWithTwoWild = HlaComplexCoverage.create(Lists.newArrayList(a1, a1, b1, b1Copy, c1, c1Copy));

        assertEquals(hom, victim.candidateRanking(Lists.newArrayList(het, hom)).get(0));
        assertEquals(het, victim.candidateRanking(Lists.newArrayList(het, homWithOneWild)).get(0));
        assertEquals(homWithOneWild, victim.candidateRanking(Lists.newArrayList(homWithOneWild, hetWithOneWild)).get(0));
        assertEquals(hetWithOneWild, victim.candidateRanking(Lists.newArrayList(homWithTwoWild, hetWithOneWild)).get(0));
    }

    */

}
