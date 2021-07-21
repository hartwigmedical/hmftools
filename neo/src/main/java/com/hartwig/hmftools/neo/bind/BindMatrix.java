package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.MAX_PEPTIDE_POSITIONS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.stats.FisherExactTest;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class BindMatrix
{
    private final String mId;
    private final int mPeptideCount;
    private final int mAminoAcidCount;

    private final double[][] mBindTotals;
    private final int[][] mBindCounts;
    private final int[][] mPredictedBinds;

    private final Map<String,double[]> mComboData;
    private int mTotal;
    private int mTotalBinds;

    private final Map<Character,Integer> mAminoAcidIndices;

    public BindMatrix(final String id, final int maxPeptidePositions)
    {
        mId = id;
        mPeptideCount = maxPeptidePositions;
        mAminoAcidCount = AMINO_ACIDS.size();
        mAminoAcidIndices = Maps.newHashMap();

        for(int i = 0; i < mAminoAcidCount; ++i)
        {
            mAminoAcidIndices.put(AMINO_ACIDS.get(i), i);
        }

        mBindTotals = new double[mAminoAcidCount][mPeptideCount];
        mBindCounts = new int[mAminoAcidCount][mPeptideCount];
        mPredictedBinds = new int[mAminoAcidCount][mPeptideCount];

        mComboData = Maps.newHashMap();
        mTotal = 0;
        mTotalBinds = 0;
    }

    public void clear()
    {
        for(int aa = 0; aa < mAminoAcidCount; ++aa)
        {
            for(int p = 0; p < mPeptideCount; ++p)
            {
                mBindCounts[aa][p] = 0;
                mBindTotals[aa][p] = 0;
            }
        }
    }

    public void logStats()
    {
        NE_LOGGER.info("id({}) total({}) binds({}) pairs({})", mId, mTotal, mTotalBinds, mComboData.size());
    }

    public void processBindData(final BindData bindData, double levelScore, boolean bindPrediction)
    {
        if(bindData.Peptide.length() != MAX_PEPTIDE_POSITIONS)
            return;

        ++mTotal;

        if(bindPrediction)
            ++mTotalBinds;

        for(int i = 0; i < bindData.Peptide.length(); ++i)
        {
            char aminoAcid = bindData.Peptide.charAt(i);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex < 0)
                continue;

            mBindTotals[aaIndex][i] += levelScore;
            ++mBindCounts[aaIndex][i];

            if(bindPrediction)
                ++mPredictedBinds[aaIndex][i];

            if(i < bindData.Peptide.length() - 1)
            {
                for(int j = i + 1; j < bindData.Peptide.length(); ++j)
                {
                    char aminoAcid2 = bindData.Peptide.charAt(j);
                    if(aminoAcidIndex(aminoAcid2) < 0)
                        continue;

                    updatePairData(aminoAcid, i, aminoAcid2, j, levelScore, bindPrediction);
                }
            }
        }
    }

    private int aminoAcidIndex(final Character aminoAcid)
    {
        Integer index = mAminoAcidIndices.get(aminoAcid);
        return index != null ? index : -1;
    }

    private static final int COMBO_DATA_COUNT = 0;
    private static final int COMBO_DATA_TOTAL = 1;
    private static final int COMBO_DATA_PRED_BINDS = 2;

    private void updatePairData(char aa1, int pos1, char aa2, int pos2, double levelScore, boolean bindPrediction)
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

        if(bindPrediction)
            ++data[COMBO_DATA_PRED_BINDS];
    }


    public static BufferedWriter initFrequencyWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,Type,AminoAcid,Peptide,Count,Score,TotalScore,PredBinds");
            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise frequency data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writeFrequencyData(final String allele, final BufferedWriter writer, final BufferedWriter probWriter)
    {
        NE_LOGGER.info("writing allele({}) frequency", allele);

        try
        {
            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                for(int p = 0; p < mPeptideCount; ++p)
                {
                    int freq = mBindCounts[aa][p];
                    int predBinds = mPredictedBinds[aa][p];
                    double totalScore = mBindTotals[aa][p];
                    double avgScore = freq > 0 ? totalScore / freq : 0;

                    writer.write(String.format("%s,Single,%s,%d,%d,%.4f,%.4f,%d",
                            allele, aminoAcid, p, freq, avgScore, totalScore, predBinds));

                    writer.newLine();
                }
            }

            // write pair-frequency data in the same format
            for(Map.Entry<String,double[]> entry : mComboData.entrySet())
            {
                String pairKey = entry.getKey();
                char aa1 = pairKey.charAt(0);
                int pos1 = Integer.parseInt(pairKey.substring(1, 2));
                char aa2 = pairKey.charAt(3);
                int pos2 = Integer.parseInt(pairKey.substring(4, 5));
                final double[] data = entry.getValue();

                int freq = (int)data[COMBO_DATA_COUNT];
                int predBinds = (int)data[COMBO_DATA_PRED_BINDS];
                double totalScore = data[COMBO_DATA_TOTAL];
                double avgScore = freq > 0 ? totalScore / freq : 0;

                writer.write(String.format("%s,Pair,%c_%c,%d_%d,%d,%.4f,%.4f,%d",
                        allele, aa1, aa2, pos1, pos2, freq, avgScore, totalScore, predBinds));

                writer.newLine();

                checkElevatedPairs(probWriter, pairKey, aa1, pos1, aa2, pos2, freq, predBinds);

            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write frequency data: {}", e.toString());
        }
    }

    public static BufferedWriter initCombProbabilityWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,Pair,PairCount,BindCount,TotalBinds,Expected,Prob,Aa1Binds,Aa1Total,Aa2Binds,Aa2Total");
            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise frequency data file({}): {}", filename, e.toString());
            return null;
        }
    }

    private void checkElevatedPairs(
            final BufferedWriter probWriter, final String pairKey, char aa1, int pos1, char aa2, int pos2,
            int pairCount, int bindingCount)
    {
        if(bindingCount == 0)
            return;

        int aa1Index = aminoAcidIndex(aa1);
        int aa1Total = mBindCounts[aa1Index][pos1];
        int aa1PredBinds = mPredictedBinds[aa1Index][pos1];

        int aa2Index = aminoAcidIndex(aa2);
        int aa2Total = mBindCounts[aa2Index][pos2];
        int aa2PredBinds = mPredictedBinds[aa2Index][pos2];

        double a1Rate = aa1Total > 0 ? aa1PredBinds / (double)mTotalBinds : 0; // was using a denom of aa1Total
        double a2Rate = aa2Total > 0 ? aa2PredBinds / (double)mTotalBinds : 0;
        double expectedPairBind = max(mTotalBinds * a1Rate * a2Rate, 0.00001);

        PoissonDistribution poisson = new PoissonDistribution(expectedPairBind);

        double poisProb = bindingCount > expectedPairBind ?
                1 - poisson.cumulativeProbability(bindingCount - 1) : poisson.cumulativeProbability(bindingCount);

        try
        {
            probWriter.write(String.format("%s,%s,%d,%d,%d,%.2f,%6.3e,%d,%d,%d,%d",
                    mId, pairKey, pairCount, bindingCount, mTotalBinds, expectedPairBind, poisProb,
                    aa1PredBinds, aa1Total, aa2PredBinds, aa2Total));

            probWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write prob data: {}", e.toString());
        }
    }
}
