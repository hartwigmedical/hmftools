package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;

import java.util.Map;

import com.google.common.collect.Maps;

public class BindScoreMatrix
{
    public final String Allele;
    public final int PeptideCount;

    private final Map<Character,Integer> mAminoAcidIndices;
    private final double[][] mBindScores; // by amino acid and position

    private static final double INVALID_SCORE = -1000;

    public BindScoreMatrix(final String allele, final int peptideLength)
    {
        Allele = allele;
        PeptideCount = peptideLength;

        int aminoAcidCount = AMINO_ACIDS.size();
        mAminoAcidIndices = Maps.newHashMap();

        for(int i = 0; i < aminoAcidCount; ++i)
        {
            mAminoAcidIndices.put(AMINO_ACIDS.get(i), i);
        }

        mBindScores = new double[aminoAcidCount][PeptideCount];
    }

    public final double[][] getBindScores() { return mBindScores; }

    public double calcScore(final String peptide)
    {
        if(peptide.length() != PeptideCount)
            return INVALID_SCORE; // for now

        double score = 0;

        for(int i = 0; i < peptide.length(); ++i)
        {
            char aminoAcid = peptide.charAt(i);

            Integer aaIndex = mAminoAcidIndices.get(aminoAcid);

            if(aaIndex == null)
                return INVALID_SCORE;

            double aaPosScore = mBindScores[aaIndex][i];

            score += aaPosScore;
        }

        return score;
    }
}
