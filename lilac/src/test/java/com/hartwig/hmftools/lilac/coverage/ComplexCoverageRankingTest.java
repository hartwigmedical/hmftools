package com.hartwig.hmftools.lilac.coverage;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.junit.Test;

public class ComplexCoverageRankingTest
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
}
