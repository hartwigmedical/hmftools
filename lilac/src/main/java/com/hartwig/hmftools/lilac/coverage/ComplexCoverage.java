package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public final class ComplexCoverage implements Comparable<ComplexCoverage>
{
    public final int TotalCoverage;
    public final int UniqueCoverage;
    public final int SharedCoverage;
    public final int WildCoverage;

    private final List<AlleleCoverage> mAlleleCoverage;

    // computed values
    private double mCohortFrequencyTotal;
    private double mScore;
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
    }

    public List<AlleleCoverage> getAlleleCoverage() { return mAlleleCoverage; }

    public List<HlaAllele> getAlleles()
    {
        return mAlleleCoverage.stream().map(x -> x.Allele).collect(Collectors.toList());
    }

    public double cohortFrequencyTotal() { return mCohortFrequencyTotal; }
    public void setCohortFrequencyTotal(double total) { mCohortFrequencyTotal = total; }

    public int homozygousCount() { return mHomozygousCount; }

    public int recoveredCount() { return mRecoveredCount; }
    public void setRecoveredCount(int count) { mRecoveredCount = count; }

    public int wildcardCount() { return mWildcardCount; }
    public void setWildcardCount(int count) { mWildcardCount = count; }

    private int calcHomozygousCount()
    {
        return EXPECTED_ALLELE_COUNT - getAlleles().size();
        /*  no repeats so just take the size relative to expected
        int hgCount = 0;
        final List<HlaAllele> alleles = getAlleles();

        for(int i = 0; i < alleles.size() - 1; i = i + 2)
        {
            if(alleles.get(i).matches(alleles.get(i + 1)))
                ++hgCount;
        }

        return hgCount;
        */
    }

    public boolean isHomozygous(final HlaAllele allele)
    {
        return getAlleles().stream().filter(x -> x.Gene.equals(allele.Gene)).count() == 1;
    }

    public void setScore(double score) { mScore = score; }
    public double getScore() { return mScore; }

    public void expandToSixAlleles()
    {
        if(mAlleleCoverage.size() == EXPECTED_ALLELE_COUNT)
            return;

        // split homozygous allele coverage, and fill any missing allele if there was zero support

        List<AlleleCoverage> existingCoverage = mAlleleCoverage.stream().collect(Collectors.toList());
        mAlleleCoverage.clear();

        for(String gene : GENE_IDS)
        {
            List<AlleleCoverage> geneCoverage = existingCoverage.stream().filter(x -> x.Allele.Gene.equals(gene)).collect(Collectors.toList());

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
        if (coverage.size() != 1)
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
                single.WildCoverage - single.WildCoverage);

        List<AlleleCoverage> newCoverage = Lists.newArrayList(first, remainder);
        return newCoverage;
    }

    public void populateMissingCoverage(final List<HlaAllele> alleles)
    {
        if(mAlleleCoverage.size() == EXPECTED_ALLELE_COUNT)
            return;

        // fill any missing allele if there was zero support
        List<AlleleCoverage> existingCoverage = mAlleleCoverage.stream().collect(Collectors.toList());
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

    public static ComplexCoverage create(final List<AlleleCoverage> alleles)
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

        final List<AlleleCoverage> sortedAlleles = alleles.stream().collect(Collectors.toList());
        Collections.sort(sortedAlleles, new AlleleCoverage.AlleleSorter());

        return new ComplexCoverage(unique, (int)round(shared), (int)round(wild), sortedAlleles);
    }

    public String toString()
    {
        return String.format("alleles(%s) coverage(%d) score(%.2f)", HlaAllele.toString(getAlleles()), TotalCoverage, mScore);
    }
}
