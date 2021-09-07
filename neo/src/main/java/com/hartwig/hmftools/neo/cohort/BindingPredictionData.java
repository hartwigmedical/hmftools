package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.gene;

import java.util.List;

public class BindingPredictionData
{
    public final String Allele;
    public final int NeId;
    public final String Peptide;

    // MhcFlurry predictions
    private double mAffinity;
    private double mAffinityPercentile;
    private double mPresentation;
    private double mPresentationPercentile;

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

        mAffinity = 0;
        mAffinityPercentile = 0;
        mPresentation = 0;
        mPresentationPercentile = 0;

        mScore = 0;
        mRankPercentile = 0;
        mLikelihood = 0;
    }

    public void setMcfPredictions(double affinity, double affinityPerc, double presentation, double presentationPerc)
    {
        mAffinity = affinity;
        mAffinityPercentile = affinityPerc;
        mPresentation = presentation;
        mPresentationPercentile = presentationPerc;
    }

    public double affinity() { return mAffinity; }
    public double presentation() { return mPresentation; }
    public double affinityPerc() { return mAffinityPercentile; }
    public double presentationPerc() { return mPresentationPercentile; }

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
        return String.format("allele(%s) neId(%d) pep(%s) affinity(%.2f) pres(%.4f)",
                Allele, NeId, Peptide, mAffinity, mPresentation);
    }

    public static BindingPredictionData fromMcfCsv(
            final String data, int alleleIndex, int neIdIndex, int peptideIndex,
            int affinityIndex, int affinityPercIndex, int presentationIndex, int presPercIndex)
    {
        final String[] items = data.split(DELIM, -1);

        String allele = items[alleleIndex].replaceAll("HLA-", "");

        BindingPredictionData bindData = new BindingPredictionData(
                allele, Integer.parseInt(items[neIdIndex]), items[peptideIndex]);

        bindData.setMcfPredictions(
                Double.parseDouble(items[affinityIndex]), Double.parseDouble(items[affinityPercIndex]) * 0.01,
                Double.parseDouble(items[presentationIndex]), Double.parseDouble(items[presPercIndex]) * 0.01);

        return bindData;

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
