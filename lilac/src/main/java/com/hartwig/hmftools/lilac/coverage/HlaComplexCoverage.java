package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_IDS;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public final class HlaComplexCoverage implements Comparable<HlaComplexCoverage>
{
    public final int TotalCoverage;
    public final int UniqueCoverage;
    public final int SharedCoverage;
    public final int WildCoverage;

    private final List<HlaAlleleCoverage> mAlleleCoverage;

    // computed values
    private double mCohortFrequencyTotal;
    private double mScore;
    private final int mHomozygousCount;
    private int mRecoveredCount;
    private int mWildcardCount;

    public HlaComplexCoverage(int uniqueCoverage, int sharedCoverage, int wildCoverage, final List<HlaAlleleCoverage> alleleCoverage)
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

    public List<HlaAlleleCoverage> getAlleleCoverage() { return mAlleleCoverage; }

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

        List<HlaAlleleCoverage> existingCoverage = mAlleleCoverage.stream().collect(Collectors.toList());
        mAlleleCoverage.clear();

        for(String gene : GENE_IDS)
        {
            List<HlaAlleleCoverage> geneCoverage = existingCoverage.stream().filter(x -> x.Allele.Gene.equals(gene)).collect(Collectors.toList());

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

    private static List<HlaAlleleCoverage> splitHomozygousCoverage(final List<HlaAlleleCoverage> coverage)
    {
        if (coverage.size() != 1)
            return coverage;

        HlaAlleleCoverage single = coverage.get(0);
        HlaAlleleCoverage first = new HlaAlleleCoverage(
                single.Allele,
                single.UniqueCoverage / 2,
                single.SharedCoverage / 2,
                single.WildCoverage / 2);

        HlaAlleleCoverage remainder = new HlaAlleleCoverage(single.Allele,
                single.UniqueCoverage - first.UniqueCoverage,
                single.SharedCoverage - first.SharedCoverage,
                single.WildCoverage - single.WildCoverage);

        List<HlaAlleleCoverage> newCoverage = Lists.newArrayList(first, remainder);
        return newCoverage;
    }

    public void populateMissingCoverage(final List<HlaAllele> alleles)
    {
        if(mAlleleCoverage.size() == EXPECTED_ALLELE_COUNT)
            return;

        // fill any missing allele if there was zero support
        List<HlaAlleleCoverage> existingCoverage = mAlleleCoverage.stream().collect(Collectors.toList());
        mAlleleCoverage.clear();

        for(HlaAllele allele : alleles)
        {
            HlaAlleleCoverage coverage = existingCoverage.stream().filter(x -> x.Allele.equals(allele)).findFirst().orElse(null);

            if(coverage == null)
            {
                mAlleleCoverage.add(new HlaAlleleCoverage(allele, 0, 0, 0));
            }
            else
            {
                mAlleleCoverage.add(coverage);
            }
        }
    }

    @Override
    public int compareTo(final HlaComplexCoverage other)
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

    public static HlaComplexCoverage create(final List<HlaAlleleCoverage> alleles)
    {
        int unique = 0;
        double shared = 0.0;
        double wild = 0.0;
        for(HlaAlleleCoverage coverage : alleles)
        {
            unique += coverage.UniqueCoverage;
            shared += coverage.SharedCoverage;
            wild += coverage.WildCoverage;
        }

        final List<HlaAlleleCoverage> sortedAlleles = alleles.stream().collect(Collectors.toList());
        Collections.sort(sortedAlleles, new HlaAlleleCoverage.AlleleSorter());

        return new HlaComplexCoverage(unique, (int)round(shared), (int)round(wild), sortedAlleles);
    }

    public String toString()
    {
        return String.format("alleles(%s) coverage(%d) score(%.2f)", HlaAllele.toString(getAlleles()), TotalCoverage, mScore);
    }
}
