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
        ComplexCoverage coverage = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0)));

        assertEquals(3, coverage.homozygousCount());

        coverage = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("A*01:03"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0)));

        assertEquals(1, coverage.homozygousCount());

        coverage = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*01:02"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("A*01:03"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*03:02"), 10, 0, 0)));

        assertEquals(0, coverage.homozygousCount());
    }

    @Test
    public void testCoverageRankSorting()
    {
        ComplexCoverage coverage1 = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*01:201"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*01:01"), 10, 0, 0)));

        coverage1.setScore(34);

        // middle 2 are sorted by their alleles numerically
        ComplexCoverage coverage2 = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*02:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*02:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*02:01"), 10, 0, 0)));

        coverage2.setScore(33);

        ComplexCoverage coverage3 = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*03:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*03:01"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*03:01"), 10, 0, 0)));

        coverage3.setScore(35);

        ComplexCoverage coverage4 = ComplexCoverage.create(Lists.newArrayList(
                new AlleleCoverage(HlaAllele.fromString("A*01:78"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("B*01:02"), 10, 0, 0),
                new AlleleCoverage(HlaAllele.fromString("C*01:02"), 10, 0, 0)));

        coverage4.setScore(34);

        List<ComplexCoverage> complexCoverages = Lists.newArrayList(coverage1, coverage2, coverage3, coverage4);

        Collections.sort(complexCoverages, new ComplexCoverageRanking.ComplexCoverageSorter());
        assertEquals(coverage3, complexCoverages.get(0));
        assertEquals(coverage4, complexCoverages.get(1));
        assertEquals(coverage1, complexCoverages.get(2));
        assertEquals(coverage2, complexCoverages.get(3));
    }
}
