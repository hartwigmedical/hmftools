package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.gene;

import java.util.List;

public class BindingPredictionData
{
    public final String Allele;
    public final int NeId;
    public final String Peptide;

    // internal scoring
    private double mScore;
    private double mRankPercentile;
    private double mLikelihood;

    public static final String DELIM = ",";

    public BindingPredictionData(final String allele, int neId, String peptide)
    {
        Allele = allele;
        NeId = neId;
        Peptide = peptide;

        mScore = 0;
        mRankPercentile = 0;
        mLikelihood = 0;
    }

    public void setScoreData(double score, double rank, double likelihood)
    {
        mScore = score;
        mRankPercentile = rank;
        mLikelihood = likelihood;
    }

    public double score() { return mScore; }
    public double rankPercentile() { return mRankPercentile; }
    public double likelihood() { return mLikelihood; }

    public String toString()
    {
        return String.format("allele(%s) neId(%d) pep(%s)", Allele, NeId, Peptide);
    }

    public static void expandHomozygous(final List<BindingPredictionData> predictions)
    {
        if(predictions.size() >= EXPECTED_ALLELE_COUNT)
            return;

        int index = 0;

        while(index < predictions.size())
        {
            BindingPredictionData coverage = predictions.get(index);
            BindingPredictionData nextCoverage = index < predictions.size() - 1 ? predictions.get(index + 1) : null;

            if(nextCoverage != null && gene(coverage.Allele).equals(gene(nextCoverage.Allele)))
            {
                // both alleles present
                index += 2;
                continue;
            }

            if(nextCoverage != null)
            {
                // must be A or B
                predictions.add(index, coverage); // replicate
                index += 2;
                continue;
            }
            else
            {
                predictions.add(coverage);
                break;
            }
        }
    }

}
