package com.hartwig.hmftools.lilac.coverage;

import static com.hartwig.hmftools.lilac.coverage.ComplexCoverageRanking.solutionComplexity;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.apache.commons.lang3.tuple.Pair;
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

    @Test
    public void testSolutionComplexitySingle()
    {
        HlaAllele allele = HlaAllele.fromString("A*01:02");
        List<String> seq = List.of("A", "A", "A", "A", "A", "A");
        Map<HlaAllele, HlaSequenceLoci> aminoAcidSequenceLookup = Stream.of(Pair.of(allele, new HlaSequenceLoci(allele, seq)))
                .collect(Collectors.toMap(Pair::getKey, Pair::getValue));

        Map<String, List<Integer>> geneExonBoundaries = Maps.newHashMap();
        geneExonBoundaries.put("HLA-A", List.of(1, 2, 100));

        List<AlleleCoverage> alleleCoverages = List.of(new AlleleCoverage(allele, 0, 0, 0));
        ComplexCoverage complexCoverage = ComplexCoverage.create(alleleCoverages);

        int actualComplexity = solutionComplexity(aminoAcidSequenceLookup, geneExonBoundaries, complexCoverage);
        int expectedComplexity = 3;

        assertEquals(expectedComplexity, actualComplexity);
    }

    @Test
    public void testSolutionComplexityThreeIdentical()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:02");
        HlaAllele allele2 = HlaAllele.fromString("A*01:03");
        HlaAllele allele3 = HlaAllele.fromString("A*02:02");

        List<String> seq = List.of("A", "A", "A", "A", "A", "A");
        Map<HlaAllele, HlaSequenceLoci> aminoAcidSequenceLookup = Stream.of(
                        Pair.of(allele1, new HlaSequenceLoci(allele1, seq)),
                        Pair.of(allele2, new HlaSequenceLoci(allele2, seq)),
                        Pair.of(allele3, new HlaSequenceLoci(allele3, seq)))
                .collect(Collectors.toMap(Pair::getKey, Pair::getValue));

        Map<String, List<Integer>> geneExonBoundaries = Maps.newHashMap();
        geneExonBoundaries.put("HLA-A", List.of(1, 2, 100));

        List<AlleleCoverage> alleleCoverages = Stream.of(allele1, allele2, allele3)
                .map(x -> new AlleleCoverage(x, 0, 0, 0))
                .toList();
        ComplexCoverage complexCoverage = ComplexCoverage.create(alleleCoverages);

        int actualComplexity = solutionComplexity(aminoAcidSequenceLookup, geneExonBoundaries, complexCoverage);
        int expectedComplexity = 3;

        assertEquals(expectedComplexity, actualComplexity);
    }

    @Test
    public void testSolutionComplexityWithMismatches()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:02");
        HlaAllele allele2 = HlaAllele.fromString("A*01:03");
        HlaAllele allele3 = HlaAllele.fromString("A*02:02");

        Map<HlaAllele, HlaSequenceLoci> aminoAcidSequenceLookup = Stream.of(
                        Pair.of(allele1, new HlaSequenceLoci(allele1, List.of("A", "A", "A", "A", "A", "A"))),
                        Pair.of(allele2, new HlaSequenceLoci(allele2, List.of("A", "B", "A", "A", "A", "A"))),
                        Pair.of(allele3, new HlaSequenceLoci(allele3, List.of("A", "A", "B", "A", "A", "A"))))
                .collect(Collectors.toMap(Pair::getKey, Pair::getValue));

        Map<String, List<Integer>> geneExonBoundaries = Maps.newHashMap();
        geneExonBoundaries.put("HLA-A", List.of(1, 2, 100));

        List<AlleleCoverage> alleleCoverages = Stream.of(allele1, allele2, allele3)
                .map(x -> new AlleleCoverage(x, 0, 0, 0))
                .toList();
        ComplexCoverage complexCoverage = ComplexCoverage.create(alleleCoverages);

        int actualComplexity = solutionComplexity(aminoAcidSequenceLookup, geneExonBoundaries, complexCoverage);
        int expectedComplexity = 5;

        assertEquals(expectedComplexity, actualComplexity);
    }

    @Test
    public void testSolutionComplexityWithWildcards()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:02");
        HlaAllele allele2 = HlaAllele.fromString("A*01:03");
        HlaAllele allele3 = HlaAllele.fromString("A*02:02");

        Map<HlaAllele, HlaSequenceLoci> aminoAcidSequenceLookup = Stream.of(
                        Pair.of(allele1, new HlaSequenceLoci(allele1, List.of("A", "A", "A", "A", "A", "A"))),
                        Pair.of(allele2, new HlaSequenceLoci(allele2, List.of("A", "B", "A", "A", "B", "A"))),
                        Pair.of(allele3, new HlaSequenceLoci(allele3, List.of("A", "A", "B", "A", "*", "A"))))
                .collect(Collectors.toMap(Pair::getKey, Pair::getValue));

        Map<String, List<Integer>> geneExonBoundaries = Maps.newHashMap();
        geneExonBoundaries.put("HLA-A", List.of(1, 2, 100));

        List<AlleleCoverage> alleleCoverages = Stream.of(allele1, allele2, allele3)
                .map(x -> new AlleleCoverage(x, 0, 0, 0))
                .toList();
        ComplexCoverage complexCoverage = ComplexCoverage.create(alleleCoverages);

        int actualComplexity = solutionComplexity(aminoAcidSequenceLookup, geneExonBoundaries, complexCoverage);
        int expectedComplexity = 5;

        assertEquals(expectedComplexity, actualComplexity);
    }

    @Test
    public void testSolutionComplexityNonUniformLengths()
    {
        HlaAllele allele1 = HlaAllele.fromString("A*01:02");
        HlaAllele allele2 = HlaAllele.fromString("A*01:03");
        HlaAllele allele3 = HlaAllele.fromString("A*02:02");

        Map<HlaAllele, HlaSequenceLoci> aminoAcidSequenceLookup = Stream.of(
                        Pair.of(allele1, new HlaSequenceLoci(allele1, List.of("A", "A", "A", "A", "A", "A"))),
                        Pair.of(allele2, new HlaSequenceLoci(allele2, List.of("A", "B", "A", "A", "B", "A"))),
                        Pair.of(allele3, new HlaSequenceLoci(allele3, List.of("A", "A", "B", "A", "*", "A", "C"))))
                .collect(Collectors.toMap(Pair::getKey, Pair::getValue));

        Map<String, List<Integer>> geneExonBoundaries = Maps.newHashMap();
        geneExonBoundaries.put("HLA-A", List.of(1, 2, 100));

        List<AlleleCoverage> alleleCoverages = Stream.of(allele1, allele2, allele3)
                .map(x -> new AlleleCoverage(x, 0, 0, 0))
                .toList();
        ComplexCoverage complexCoverage = ComplexCoverage.create(alleleCoverages);

        int actualComplexity = solutionComplexity(aminoAcidSequenceLookup, geneExonBoundaries, complexCoverage);
        int expectedComplexity = 7;

        assertEquals(expectedComplexity, actualComplexity);
    }
}
