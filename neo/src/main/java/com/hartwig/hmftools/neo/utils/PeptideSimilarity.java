package com.hartwig.hmftools.neo.utils;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

public class PeptideSimilarity
{
    // sourced from external tool eg Prime
    public final String Allele;
    public final String Peptide;
    public final String WildtypePeptide;

    // closest in the proteome regardless of binding
    private String mTopPeptide;
    private double mTopSimiliarity;
    private TranscriptAminoAcids mTopTransData;

    private String mNearestBinder;
    private double mNearestBinderLikelihoodRank;

    public PeptideSimilarity(final String peptide, final String allele, final String wtPeptide)
    {
        Peptide = peptide;
        Allele = allele;
        WildtypePeptide = wtPeptide;

        mTopTransData = null;
        mTopSimiliarity = -1;
        mTopPeptide = "";
        mNearestBinder = "";
        mNearestBinderLikelihoodRank = -1;
    }

    public void setTopSimilarity(final String peptide, double similiarity, final TranscriptAminoAcids transData)
    {
        mTopPeptide = peptide;
        mTopSimiliarity = similiarity;
        mTopTransData = transData;
    }

    public String topPeptide() { return mTopPeptide; }
    public double topSimiliarity() { return mTopSimiliarity; }
    public TranscriptAminoAcids topTransData() { return mTopTransData; }

    public void setTopBinder(final String peptide, double rank)
    {
        mNearestBinder = peptide;
        mNearestBinderLikelihoodRank = rank;
    }

    public String nearestBinder() { return mNearestBinder; }
    public double nearestBinderLikelihoodRank() { return mNearestBinderLikelihoodRank; }

    public String positionDiffs(final String otherPeptide) { return positionDiffs(Peptide, otherPeptide); }

    public String toString() { return String.format("allele(%s) peptide(%s)", Allele, Peptide); }

    public static String positionDiffs(final String origPeptide, final String otherPeptide)
    {
        if(origPeptide.equals(otherPeptide))
            return "";

        StringJoiner sj = new StringJoiner(";");
        for(int i = 0; i < otherPeptide.length(); ++i)
        {
            if(otherPeptide.charAt(i) != origPeptide.charAt(i))
                sj.add(String.valueOf(i));
        }

        return sj.toString();
    }
}
