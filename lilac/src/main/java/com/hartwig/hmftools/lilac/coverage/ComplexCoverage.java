package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;

public final class ComplexCoverage implements Comparable<ComplexCoverage>
{
    public final int TotalCoverage;
    public final int UniqueCoverage;
    public final int SharedCoverage;
    public final int WildCoverage;

    private final List<AlleleCoverage> mAlleleCoverage;

    // computed values
    private double mCohortFrequencyPenalty;
    private double mCohortFrequencyTotal;
    private double mScore;
    private double mComplexityPenalty;
    private int mComplexity;
    private final int mHomozygousCount;
    private int mRecoveredCount;
    private int mWildcardCount;

    public ComplexCoverage(int uniqueCoverage, int sharedCoverage, int wildCoverage, final List<AlleleCoverage> alleleCoverage)
    {
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;

        mAlleleCoverage = alleleCoverage;

        mHomozygousCount = calcHomozygousCount();
        mRecoveredCount = 0;
        mWildcardCount = 0;
        mCohortFrequencyTotal = 0;
        mScore = 0;
        mComplexityPenalty = 0;
        mComplexity = 0;
        mCohortFrequencyPenalty = 0;
    }

    public List<AlleleCoverage> getAlleleCoverage() { return mAlleleCoverage; }

    public List<HlaAllele> getAlleles()
    {
        return mAlleleCoverage.stream().map(x -> x.Allele).collect(Collectors.toList());
    }

    public HlaComplex toComplex() { return new HlaComplex(getAlleles()); }

    public double cohortFrequencyTotal() { return mCohortFrequencyTotal; }
    public void setCohortFrequencyTotal(double total) { mCohortFrequencyTotal = total; }

    public int homozygousCount() { return mHomozygousCount; }

    public int recoveredCount() { return mRecoveredCount; }
    public void setRecoveredCount(int count) { mRecoveredCount = count; }

    public int wildcardCount() { return mWildcardCount; }
    public void setWildcardCount(int count) { mWildcardCount = count; }

    private int calcHomozygousCount()
    {
        // for unit tests
        if(GENE_CACHE == null)
            return 0;

        return GENE_CACHE.ExpectAlleleCount - getAlleles().size();
    }

    public boolean isHomozygous(final HlaAllele allele)
    {
        return getAlleles().stream().filter(x -> x.Gene == allele.Gene).count() == 1;
    }

    public void setScore(double score) { mScore = score; }
    public double getScore() { return mScore; }
    public void setComplexityPenalty(double complexityPenalty) { mComplexityPenalty = complexityPenalty; }
    public double getComplexityPenalty() { return mComplexityPenalty; }
    public void setComplexity(int complexity) { mComplexity = complexity; }
    public int getComplexity() { return mComplexity; }
    public void setCohortFrequencyPenalty(double cohortFrequencyPenalty) { mCohortFrequencyPenalty = cohortFrequencyPenalty; }
    public double getCohortFrequencyPenalty() { return mCohortFrequencyPenalty; }

    public void expandToSixAlleles()
    {
        if(mAlleleCoverage.size() == GENE_CACHE.ExpectAlleleCount)
            return;

        // split homozygous allele coverage, and fill any missing allele if there was zero support

        List<AlleleCoverage> existingCoverage = mAlleleCoverage.stream().toList();
        mAlleleCoverage.clear();

        for(HlaGene gene : GENE_CACHE.GeneNames)
        {
            List<AlleleCoverage> geneCoverage = existingCoverage.stream()
                    .filter(x -> x.Allele.Gene == gene)
                    .collect(Collectors.toList());

            if(geneCoverage.size() == 2)
            {
                mAlleleCoverage.addAll(geneCoverage);
            }
            else
            {
                mAlleleCoverage.addAll(splitHomozygousCoverage(geneCoverage));
            }
        }
    }

    private static List<AlleleCoverage> splitHomozygousCoverage(final List<AlleleCoverage> coverage)
    {
        if(coverage.size() != 1)
            return coverage;

        AlleleCoverage single = coverage.get(0);
        AlleleCoverage first = new AlleleCoverage(
                single.Allele,
                single.UniqueCoverage / 2,
                single.SharedCoverage / 2,
                single.WildCoverage / 2);

        AlleleCoverage remainder = new AlleleCoverage(single.Allele,
                single.UniqueCoverage - first.UniqueCoverage,
                single.SharedCoverage - first.SharedCoverage,
                0.0);

        List<AlleleCoverage> newCoverage = Lists.newArrayList(first, remainder);
        return newCoverage;
    }

    public void populateMissingCoverage(final LinkedHashSet<HlaAllele> alleles)
    {
        if(mAlleleCoverage.size() == GENE_CACHE.ExpectAlleleCount)
            return;

        // fill any missing allele if there was zero support
        List<AlleleCoverage> existingCoverage = mAlleleCoverage.stream().toList();
        mAlleleCoverage.clear();

        for(HlaAllele allele : alleles)
        {
            AlleleCoverage coverage = existingCoverage.stream().filter(x -> x.Allele.equals(allele)).findFirst().orElse(null);

            if(coverage == null)
            {
                mAlleleCoverage.add(new AlleleCoverage(allele, 0, 0, 0));
            }
            else
            {
                mAlleleCoverage.add(coverage);
            }
        }
    }

    @Override
    public int compareTo(final ComplexCoverage other)
    {
        if(TotalCoverage != other.TotalCoverage)
            return TotalCoverage > other.TotalCoverage ? 1 : -1;

        if(SharedCoverage != other.SharedCoverage)
            return SharedCoverage > other.SharedCoverage ? 1 : -1;

        if(WildCoverage != other.WildCoverage)
            return WildCoverage > other.WildCoverage ? 1 : -1;

        if(UniqueCoverage != other.UniqueCoverage)
            return UniqueCoverage > other.UniqueCoverage ? 1 : -1;

        return 0;
    }

    public static ComplexCoverage create(final Iterable<AlleleCoverage> alleles)
    {
        int unique = 0;
        double shared = 0.0;
        double wild = 0.0;
        for(AlleleCoverage coverage : alleles)
        {
            unique += coverage.UniqueCoverage;
            shared += coverage.SharedCoverage;
            wild += coverage.WildCoverage;
        }

        final List<AlleleCoverage> sortedAlleles = Lists.newArrayList(alleles);
        Collections.sort(sortedAlleles, new AlleleCoverage.AlleleSorter());

        return new ComplexCoverage(unique, (int) round(shared), (int) round(wild), sortedAlleles);
    }

    @Override
    public String toString()
    {
        return String.format("alleles(%s) coverage(%d) score(%.2f) complexityPenalty(%.2f)", HlaAllele.toString(getAlleles()), TotalCoverage, mScore, mComplexityPenalty);
    }
}
