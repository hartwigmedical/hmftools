package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.MAX_PEPTIDE_POSITIONS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;

public class BindMatrix
{
    private final int mPeptideCount;
    private final int mAminoAcidCount;

    private final double[][] mBindingTotals;
    private final int[][] mBindingCounts;

    private final Map<Character,Integer> mAminoAcidIndices;

    public BindMatrix(final int maxPeptidePositions)
    {
        mPeptideCount = maxPeptidePositions;
        mAminoAcidCount = AMINO_ACIDS.size();
        mAminoAcidIndices = Maps.newHashMap();

        for(int i = 0; i < mAminoAcidCount; ++i)
        {
            mAminoAcidIndices.put(AMINO_ACIDS.get(i), i);
        }

        mBindingTotals = new double[mAminoAcidCount][mPeptideCount];
        mBindingCounts = new int[mAminoAcidCount][mPeptideCount];
    }

    public void clear()
    {
        for(int aa = 0; aa < mAminoAcidCount; ++aa)
        {
            for(int p = 0; p < mPeptideCount; ++p)
            {
                mBindingCounts[aa][p] = 0;
                mBindingTotals[aa][p] = 0;
            }
        }
    }

    public void processBindData(final BindData bindData, double levelScore)
    {
        if(bindData.Peptide.length() != MAX_PEPTIDE_POSITIONS)
            return;

        for(int i = 0; i < bindData.Peptide.length(); ++i)
        {
            char aminoAcid = bindData.Peptide.charAt(i);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex < 0)
                continue;

            mBindingTotals[aaIndex][i] += levelScore;
            ++mBindingCounts[aaIndex][i];
        }
    }

    private int aminoAcidIndex(final Character aminoAcid)
    {
        Integer index = mAminoAcidIndices.get(aminoAcid);
        return index != null ? index : -1;
    }

    public static BufferedWriter initialiseFrequencyWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,AminoAcid,Peptide,Count,Score,TotalScore");
            writer.newLine();

            return writer;
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to initialise frequency data file({}): {}", filename, e.toString());
            return null;
        }
    }

    public void writeFrequencyData(final String allele, final BufferedWriter writer)
    {
        NE_LOGGER.info("writing allele({}) frequency", allele);

        try
        {
            for(int aa = 0; aa < mAminoAcidCount; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                for(int p = 0; p < mPeptideCount; ++p)
                {
                    int freq = mBindingCounts[aa][p];
                    double totalScore = mBindingTotals[aa][p];
                    double avgScore = freq > 0 ? totalScore / freq : 0;

                    writer.write(String.format("%s,%s,%s,%d,%.4f,%.4f",
                            allele, aminoAcid, p + 1, freq, avgScore, totalScore));

                    writer.newLine();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write frequency data: {}", e.toString());
        }
    }

    public void writeMatrixResults(final String outputDir, final String filename)
    {

    }

}
