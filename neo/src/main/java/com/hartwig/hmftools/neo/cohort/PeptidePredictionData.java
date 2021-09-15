package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.bind.BindScorer.INVALID_CALC;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.gene;

import java.util.List;

public class PeptidePredictionData
{
    public final int NeId;
    public final String Allele;
    public final String Peptide;

    // internal scoring
    private double mScore;
    private double mRankPercentile;
    private double mLikelihood;
    private double mLikelihoodRank;

    public static final String DELIM = ",";

    public PeptidePredictionData(int neId, final String allele, String peptide)
    {
        Allele = allele;
        NeId = neId;
        Peptide = peptide;

        mScore = INVALID_CALC;
        mRankPercentile = INVALID_CALC;
        mLikelihood = INVALID_CALC;
        mLikelihoodRank = INVALID_CALC;
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

    public static void expandHomozygous(final List<PeptidePredictionData> predictions)
    {
        if(predictions.size() >= EXPECTED_ALLELE_COUNT)
            return;

        int index = 0;

        while(index < predictions.size())
        {
            PeptidePredictionData coverage = predictions.get(index);
            PeptidePredictionData nextCoverage = index < predictions.size() - 1 ? predictions.get(index + 1) : null;

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
