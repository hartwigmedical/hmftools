package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class ScoringData
{
    private final String mAllele;
    private final CalcConstants mCalcConstants;
    private final int mPeptideCount;
    private final int mAminoAcidCount;

    private final int[][] mObservations; // data points with a given amino acid and peptide
    private final double[][] mBindScoreTotals; // log-score total
    private final double[][] mActualBinds; // observations with affinity below configured binding threshold

    private final Map<String,double[]> mComboData;
    private int mTotal;
    private int mTotalBinds; // vs high affinity threshold
    private double mCalcTotalBinds; // vs high affinity threshold

    public ScoringData(final String allele, final int peptideLength, final CalcConstants calcConstants)
    {
        mAllele = allele;
        mCalcConstants = calcConstants;
        mPeptideCount = peptideLength;
        mAminoAcidCount = AMINO_ACIDS.size();

        mObservations = new int[mAminoAcidCount][mPeptideCount];
        mBindScoreTotals = new double[mAminoAcidCount][mPeptideCount];
        mActualBinds = new double[mAminoAcidCount][mPeptideCount];

        mComboData = Maps.newHashMap();
        mTotal = 0;
        mTotalBinds = 0;
        mCalcTotalBinds = 0;
    }

    public void logStats()
    {
        NE_LOGGER.info("allele({}) peptideLen({}) total({}) binds({}) pairs({})",
                mAllele, mPeptideCount, mTotal, mTotalBinds, mComboData.size());
    }

    public void processBindData(final BindData bindData, boolean calcPairs)
    {
        if(bindData.Peptide.length() != mPeptideCount || !bindData.Allele.equals(mAllele))
            return;

        double levelScore = mCalcConstants.deriveLevelScore(bindData.Affinity);
        double bindPerc = mCalcConstants.deriveAffinityPercent(bindData.Affinity);

        ++mTotal;
        mCalcTotalBinds += bindPerc;

        boolean actualBind = bindData.Affinity < mCalcConstants.BindingAffinityHigh;
        boolean predictedBind = bindData.PredictedAffinity < mCalcConstants.BindingAffinityHigh && actualBind;

        if(actualBind)
            ++mTotalBinds;

        for(int pos = 0; pos < bindData.Peptide.length(); ++pos)
        {
            char aminoAcid = bindData.Peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            ++mObservations[aaIndex][pos];
            mBindScoreTotals[aaIndex][pos] += levelScore;
            mActualBinds[aaIndex][pos] += bindPerc;

            if(calcPairs && bindData.isTraining())
            {
                if(pos < bindData.Peptide.length() - 1)
                {
                    for(int pos2 = pos + 1; pos2 < bindData.Peptide.length(); ++pos2)
                    {
                        char aminoAcid2 = bindData.Peptide.charAt(pos2);
                        if(aminoAcidIndex(aminoAcid2) == INVALID_AMINO_ACID)
                            continue;

                        updatePairData(aminoAcid, pos, aminoAcid2, pos2, levelScore, actualBind, predictedBind);
                    }
                }
            }
        }
    }

    public BindScoreMatrix createMatrix(final AminoAcidFrequency aminoAcidFrequency)
    {
        NE_LOGGER.debug("creating allele({}) peptideLength({}) matrix data", mAllele, mPeptideCount);

        BindScoreMatrix matrix = new BindScoreMatrix(mAllele, mPeptideCount);
        final double[][] data = matrix.getBindScores();

        for(int aa = 0; aa < mAminoAcidCount; ++aa)
        {
            char aminoAcid = AMINO_ACIDS.get(aa);
            double aaFrequency = aminoAcidFrequency.getAminoAcidFrequency(aminoAcid);

            for(int pos = 0; pos < mPeptideCount; ++pos)
            {
                double totalBinds = mActualBinds[aa][pos];
                double freqPerc = totalBinds / mCalcTotalBinds;
                freqPerc = max(freqPerc, MIN_OBSERVED_AA_POS_FREQ);

                // Score = log(max(FreqPerc,0.001) / AminoAcidFreq,2)
                double score = log(2, freqPerc / aaFrequency);

                data[aa][pos] = score;
            }
        }

        return matrix;
    }

    public static BufferedWriter initMatrixWriter(final String filename, int peptideLength)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,AminoAcid");

            for(int i = 0; i < peptideLength; ++i)
            {
                writer.write(String.format(",P%d", i));
            }

            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise matrix data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writeMatrixData(final BufferedWriter writer, final BindScoreMatrix matrix, int maxPeptideLength)
    {
        NE_LOGGER.debug("writing allele({}) matrix data", mAllele);

        final double[][] data = matrix.getBindScores();

        try
        {
            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                writer.write(String.format("%s,%d,%c", mAllele, mPeptideCount, aminoAcid));

                for(int pos = 0; pos < maxPeptideLength; ++pos)
                {
                    if(pos >= mPeptideCount)
                    {
                        writer.write(",0.0");
                        continue;
                    }

                    writer.write(String.format(",%.4f", data[aa][pos]));
                }

                writer.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write matrix data: {}", e.toString());
        }
    }

    public static BufferedWriter initFrequencyWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,AminoAcid,PeptidePos,Count,Score,TotalScore,ActualBinds");
            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise frequency data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writeFrequencyData(final BufferedWriter writer)
    {
        NE_LOGGER.debug("writing allele({}) frequency", mAllele);

        try
        {
            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                for(int p = 0; p < mPeptideCount; ++p)
                {
                    int freq = mObservations[aa][p];
                    double totalBinds = mActualBinds[aa][p];
                    double totalScore = mBindScoreTotals[aa][p];
                    double avgScore = freq > 0 ? totalScore / freq : 0;

                    writer.write(String.format("%s,%d,%s,%d,%d,%.4f,%.4f,%.2f",
                            mAllele, mPeptideCount, aminoAcid, p, freq, avgScore, totalScore, totalBinds));

                    writer.newLine();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write frequency data: {}", e.toString());
        }
    }

    private static final int COMBO_DATA_COUNT = 0;
    private static final int COMBO_DATA_TOTAL = 1;
    private static final int COMBO_DATA_ACTUAL_BINDS = 2;
    private static final int COMBO_DATA_PRED_BINDS = 3;

    private void updatePairData(
            char aa1, int pos1, char aa2, int pos2, double levelScore, boolean actualBind, boolean predictedBind)
    {
        String pairKey = String.format("%c%d_%c%d", aa1, pos1, aa2, pos2);
        double[] data = mComboData.get(pairKey);

        if(data == null)
        {
            data = new double[COMBO_DATA_PRED_BINDS + 1];
            mComboData.put(pairKey, data);
        }

        ++data[COMBO_DATA_COUNT];
        data[COMBO_DATA_TOTAL] += levelScore;

        if(actualBind)
            ++data[COMBO_DATA_ACTUAL_BINDS];

        if(predictedBind)
            ++data[COMBO_DATA_PRED_BINDS];
    }

    public static BufferedWriter initPairDataWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,AminoAcids,Positions,Count,Score,TotalScore,ActualBinds,PredictedBinds");
            writer.write(",TotalBinds,Expected,Prob,Aa1Binds,Aa1Total,Aa2Binds,Aa2Total");
            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise pair data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writePairData(final String allele, final BufferedWriter writer)
    {
        NE_LOGGER.info("writing allele({}) pair probabilities", allele);

        try
        {
            for(Map.Entry<String,double[]> entry : mComboData.entrySet())
            {
                String pairKey = entry.getKey();
                char aa1 = pairKey.charAt(0);
                int pos1 = Integer.parseInt(pairKey.substring(1, 2));
                char aa2 = pairKey.charAt(3);
                int pos2 = Integer.parseInt(pairKey.substring(4, 5));
                final double[] data = entry.getValue();

                int pairCount = (int)data[COMBO_DATA_COUNT];
                int predictedBinds = (int)data[COMBO_DATA_PRED_BINDS];
                int actualBinds = (int)data[COMBO_DATA_ACTUAL_BINDS];
                double totalScore = data[COMBO_DATA_TOTAL];
                double avgScore = pairCount > 0 ? totalScore / pairCount : 0;

                int aa1Index = aminoAcidIndex(aa1);
                int aa1Total = mObservations[aa1Index][pos1];
                double aa1ActBinds = mActualBinds[aa1Index][pos1];

                int aa2Index = aminoAcidIndex(aa2);
                int aa2Total = mObservations[aa2Index][pos2];
                double aa2ActBinds = mActualBinds[aa2Index][pos2];

                double a1Rate = aa1Total > 0 ? aa1ActBinds / (double)mTotalBinds : 0; // was using a denom of aa1Total
                double a2Rate = aa2Total > 0 ? aa2ActBinds / (double)mTotalBinds : 0;
                double expectedPairBind = max(mTotalBinds * a1Rate * a2Rate, 0.00001);

                PoissonDistribution poisson = new PoissonDistribution(expectedPairBind);

                double poisProb = actualBinds > expectedPairBind ?
                        1 - poisson.cumulativeProbability(actualBinds - 1) : poisson.cumulativeProbability(actualBinds);

                writer.write(String.format("%s,%c_%c,%d_%d,%d,%.4f,%.4f,%d,%d",
                        mAllele, aa1, aa2, pos1, pos2, pairCount, avgScore, totalScore, actualBinds, predictedBinds));

                writer.write(String.format(",%d,%.2f,%6.3e,%.2f,%d,%.2f,%d",
                        mTotalBinds, expectedPairBind, poisProb, aa1ActBinds, aa1Total, aa2ActBinds, aa2Total));

                writer.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write pair frequency & probability data: {}", e.toString());
        }
    }
}
