package com.hartwig.hmftools.lilac.coverage;

import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

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

    public HlaComplexCoverage(int uniqueCoverage, int sharedCoverage, int wildCoverage, final List<HlaAlleleCoverage> alleleCoverage)
    {
        UniqueCoverage = uniqueCoverage;
        SharedCoverage = sharedCoverage;
        WildCoverage = wildCoverage;
        TotalCoverage = UniqueCoverage + SharedCoverage + WildCoverage;

        mAlleleCoverage = alleleCoverage;

        mHomozygousCount = calcHomozygousCount();
        mRecoveredCount = 0;
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

    public void setScore(double score) { mScore = score; }
    public double getScore() { return mScore; }

    public final HlaComplexCoverage expandToSixAlleles()
    {
        return create(HlaAlleleCoverage.expand(mAlleleCoverage));
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

    public static class TotalCoverageSorter implements Comparator<HlaComplexCoverage>
    {
        // sorts by total coverage descending
        public int compare(final HlaComplexCoverage first, final HlaComplexCoverage second)
        {
            if(first.TotalCoverage != second.TotalCoverage)
                return first.TotalCoverage < second.TotalCoverage ? 1 : -1;

            return 0;
        }
    }

}
