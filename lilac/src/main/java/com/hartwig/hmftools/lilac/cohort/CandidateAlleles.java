package com.hartwig.hmftools.lilac.cohort;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.lilac.hla.HlaAllele;

public class CandidateAlleles
{
    public final int Id;
    public final int TotalCoverage;
    public final List<HlaAllele> Alleles;

    private double mCohortFrequencyTotal;
    private final int mHomozygousCount;

    private static int TOTAL_COVERAGE_DENOM = 1000;

    public CandidateAlleles(final int candidateId, final List<HlaAllele> alleles, int totalCoverage)
    {
        Id = candidateId;
        TotalCoverage = totalCoverage;
        Alleles = alleles;

        mHomozygousCount = homozygousCount();
    }

    public void setCohortFrequencyTotal(double total) { mCohortFrequencyTotal = total; }

    public int homozygousCount()
    {
        int hgCount = 0;

        for(int i = 0; i < 5; i = i + 2)
        {
            if(Alleles.get(i).matches(Alleles.get(i + 1)))
                ++hgCount;
        }

        return hgCount;
    }

    public double calcScore()
    {
        return TotalCoverage
                + mCohortFrequencyTotal * 1.5 * TotalCoverage / TOTAL_COVERAGE_DENOM
                + mHomozygousCount * 4.5 * TotalCoverage / (double)TOTAL_COVERAGE_DENOM;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(",");
        sj.add("CandidateId");
        sj.add("TotalCoverage");
        sj.add("FrequencyCalc");
        sj.add("HomozygousCount");
        sj.add("Score");
        sj.add("A1");
        sj.add("A2");
        sj.add("B1");
        sj.add("B2");
        sj.add("C1");
        sj.add("C2");
        sj.add("IsTruthSet");
        sj.add("DiffFromTop");
        return sj.toString();
    }

    public String toCsv(final double topScore, final List<HlaAllele> truthSetAlleles)
    {
        StringJoiner sj = new StringJoiner(",");
        sj.add(String.valueOf(Id));
        sj.add(String.valueOf(TotalCoverage));
        sj.add(String.format("%.4f", mCohortFrequencyTotal));
        sj.add(String.format("%d", mHomozygousCount));

        double score = calcScore();
        double diffFromTop = topScore - score;

        sj.add(String.format("%.1f", calcScore()));

        boolean isTruthSet = true;

        for(int i = 0; i < 6; ++i)
        {
            HlaAllele allele = Alleles.get(i);
            sj.add(allele.toString());

            if(!allele.matches(truthSetAlleles.get(i)))
                isTruthSet = false;
        }

        sj.add(String.valueOf(isTruthSet));
        sj.add(String.format("%.4f", diffFromTop));

        return sj.toString();
    }

}
