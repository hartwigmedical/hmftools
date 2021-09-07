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
    public static void updatePairData(final Map<String,Integer> comboDataMap, char aa1, int pos1, char aa2, int pos2)
    {
        String pairKey = String.format("%c%d_%c%d", aa1, pos1, aa2, pos2);
        Integer count = comboDataMap.get(pairKey);

        if(count == null)
            comboDataMap.put(pairKey,  1);
        else
            comboDataMap.put(pairKey,  count + 1);
    }

    public static BufferedWriter initPairDataWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,AminoAcids,Positions,PairBindCount");
            writer.write(",TotalBinds,Expected,Prob,Aa1Binds,Aa2Binds");
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
        final Map<String,Integer> comboData = bindCountData.getComboData();
        int totalBinds = bindCountData.totalBindCount();

        try
        {
            for(Map.Entry<String,Integer> entry : comboData.entrySet())
            {
                String pairKey = entry.getKey();
                char aa1 = pairKey.charAt(0);
                int pos1 = Integer.parseInt(pairKey.substring(1, 2));
                char aa2 = pairKey.charAt(3);
                int pos2 = Integer.parseInt(pairKey.substring(4, 5));
                int pairBindCount = entry.getValue();

                int aa1Index = aminoAcidIndex(aa1);
                double aa1ActBinds = bindCounts[aa1Index][pos1];

                int aa2Index = aminoAcidIndex(aa2);
                double aa2ActBinds = bindCounts[aa2Index][pos2];

                double a1Rate = aa1ActBinds / (double)totalBinds;
                double a2Rate = aa2ActBinds / (double)totalBinds;
                double expectedPairBind = max(totalBinds * a1Rate * a2Rate, 0.00001);

                PoissonDistribution poisson = new PoissonDistribution(expectedPairBind);

                double poisProb = pairBindCount > expectedPairBind ?
                        1 - poisson.cumulativeProbability(pairBindCount - 1) : poisson.cumulativeProbability(pairBindCount);

                writer.write(String.format("%s,%c_%c,%d_%d,%d,%d,%d",
                        bindCountData.Allele, aa1, aa2, pos1, pos2, pairBindCount));

                writer.write(String.format(",%d,%.2f,%6.3e,%d,%d",
                        totalBinds, expectedPairBind, poisProb, aa1ActBinds, aa2ActBinds));

                writer.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write pair frequency & probability data: {}", e.toString());
        }
    }

}
