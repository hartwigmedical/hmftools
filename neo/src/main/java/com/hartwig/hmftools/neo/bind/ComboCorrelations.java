package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import org.apache.commons.math3.distribution.PoissonDistribution;

public final class ComboCorrelations
{
    private static final int COMBO_DATA_COUNT = 0;
    private static final int COMBO_DATA_TOTAL = 1;
    private static final int COMBO_DATA_ACTUAL_BINDS = 2;
    private static final int COMBO_DATA_PRED_BINDS = 3;

    public static void updatePairData(
            final Map<String,double[]> comboDataMap,
            char aa1, int pos1, char aa2, int pos2, double levelScore, boolean actualBind, boolean predictedBind)
    {
        String pairKey = String.format("%c%d_%c%d", aa1, pos1, aa2, pos2);
        double[] data = comboDataMap.get(pairKey);

        if(data == null)
        {
            data = new double[COMBO_DATA_PRED_BINDS + 1];
            comboDataMap.put(pairKey, data);
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

    public static void writePairData(
            final BufferedWriter writer, final BindCountData bindCountData)
    {
        NE_LOGGER.info("writing allele({}) pair probabilities", bindCountData.Allele);

        final double[][] bindCounts = bindCountData.getBindCounts();
        final int[][] observations = bindCountData.getObservations();
        final Map<String,double[]> comboData = bindCountData.getComboData();
        int totalBinds = bindCountData.totalBindCount();

        try
        {
            for(Map.Entry<String,double[]> entry : comboData.entrySet())
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
                int aa1Total = observations[aa1Index][pos1];
                double aa1ActBinds = bindCounts[aa1Index][pos1];

                int aa2Index = aminoAcidIndex(aa2);
                int aa2Total = observations[aa2Index][pos2];
                double aa2ActBinds = bindCounts[aa2Index][pos2];

                double a1Rate = aa1Total > 0 ? aa1ActBinds / (double)totalBinds : 0; // was using a denom of aa1Total
                double a2Rate = aa2Total > 0 ? aa2ActBinds / (double)totalBinds : 0;
                double expectedPairBind = max(totalBinds * a1Rate * a2Rate, 0.00001);

                PoissonDistribution poisson = new PoissonDistribution(expectedPairBind);

                double poisProb = actualBinds > expectedPairBind ?
                        1 - poisson.cumulativeProbability(actualBinds - 1) : poisson.cumulativeProbability(actualBinds);

                writer.write(String.format("%s,%c_%c,%d_%d,%d,%.4f,%.4f,%d,%d",
                        bindCountData.Allele, aa1, aa2, pos1, pos2, pairCount, avgScore, totalScore, actualBinds, predictedBinds));

                writer.write(String.format(",%d,%.2f,%6.3e,%.2f,%d,%.2f,%d",
                        bindCounts, expectedPairBind, poisProb, aa1ActBinds, aa1Total, aa2ActBinds, aa2Total));

                writer.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write pair frequency & probability data: {}", e.toString());
        }
    }

}
