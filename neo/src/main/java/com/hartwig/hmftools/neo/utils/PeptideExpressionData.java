package com.hartwig.hmftools.neo.utils;

import static com.hartwig.hmftools.neo.bind.BindCommon.ITEM_DELIM;

import java.util.List;
import java.util.StringJoiner;

class PeptideExpressionData
{
    public final String Allele;
    public final String Peptide;
    public final String Source;
    public final double LikelihoodRank;
    public final List<String> Transcripts;

    private double mTpmTotal;
    private boolean mHasTpm;

    public static final String SOURCE_PROTEOME = "Proteome";
    public static final String SOURCE_VALIDATION = "Validation";

    public PeptideExpressionData(
            final String allele, final String peptide, final String source, final double likelihoodRank, final List<String> transcripts)
    {
        Allele = allele;
        Peptide = peptide;
        Source = source;
        LikelihoodRank = likelihoodRank;
        Transcripts = transcripts;

        mTpmTotal = 0;
        mHasTpm = false;
    }

    public void addTpm(final Double tpm)
    {
        if(tpm == null)
            return;

        mHasTpm = true;
        mTpmTotal += tpm;
    }

    public boolean hasTpm() { return mHasTpm; }
    public double tpm() { return mTpmTotal; }

    public String toString() { return String.format("allele(%s) peptide(%s) source(%s)", Allele, Peptide, Source); }

    public String transcripts()
    {
        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        Transcripts.forEach(x -> sj.add(x));
        return sj.toString();
    }
}
