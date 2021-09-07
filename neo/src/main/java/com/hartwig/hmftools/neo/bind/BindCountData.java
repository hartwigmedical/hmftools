package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_ALLELE_WEIGHTED;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_BIND_COUNTS;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_LENGTH_WEIGHTED;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_NOISE;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.BindCommon.COUNT_DATA_TYPES;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.MatrixUtils;

public class BindCountData
{
    public final String Allele;
    public final int PeptideLength;

    private final double[][] mBindCounts; // binding observations
    private final double[][] mNoiseCounts; // adjusted for noise

    private final double[][] mWeightedCounts; // weighted across the peptide lengths for this allele, padded to reference length
    private final double[] mAdjustedPosTotals; // weighted totals per position from all lengths
    private final double[][] mFinalWeightedCounts; // weighted across all other alleles, padded to reference length

    private final Map<String,Integer> mComboCounts;
    private int mTotalBinds;

    public BindCountData(final String allele, final int peptideLength)
    {
        Allele = allele;
        PeptideLength = peptideLength;
        int aminoAcidCount = AMINO_ACID_COUNT;

        mBindCounts = new double[aminoAcidCount][PeptideLength];
        mNoiseCounts = new double[aminoAcidCount][PeptideLength];
        mWeightedCounts = new double[aminoAcidCount][PeptideLength];
        mFinalWeightedCounts = new double[aminoAcidCount][PeptideLength];
        mAdjustedPosTotals = new double[PeptideLength];

        mComboCounts = Maps.newHashMap();
        mTotalBinds = 0;
    }

    public final double[][] getBindCounts() { return mBindCounts; }
    public Map<String,Integer> getComboData() { return mComboCounts; }
    public final double[][] getNoiseCounts() { return mNoiseCounts; }
    public final double[][] getWeightedCounts() { return mWeightedCounts; }
    public final double[] getAdjustedPosTotals() { return mAdjustedPosTotals; }
    public final double[][] getFinalWeightedCounts() { return mFinalWeightedCounts; }
    public int totalBindCount() { return mTotalBinds; }

    public void logStats()
    {
        NE_LOGGER.debug("allele({}) peptideLen({}) totalBinds({}) pairs({})",
                Allele, PeptideLength, mTotalBinds, mComboCounts.size());
    }

    public void processBindData(final BindData bindData, boolean calcPairs)
    {
        if(bindData.Peptide.length() != PeptideLength || !bindData.Allele.equals(Allele))
            return;

        ++mTotalBinds; // every observation means the peptide does bind

        for(int pos = 0; pos < bindData.Peptide.length(); ++pos)
        {
            char aminoAcid = bindData.Peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            ++mBindCounts[aaIndex][pos];

            if(calcPairs)
            {
                if(pos < bindData.Peptide.length() - 1)
                {
                    for(int pos2 = pos + 1; pos2 < bindData.Peptide.length(); ++pos2)
                    {
                        char aminoAcid2 = bindData.Peptide.charAt(pos2);
                        if(aminoAcidIndex(aminoAcid2) == INVALID_AMINO_ACID)
                            continue;

                        ComboCorrelations.updatePairData(mComboCounts, aminoAcid, pos, aminoAcid2, pos2);
                    }
                }
            }
        }
    }

    public void processBindingPeptide(final String peptide)
    {
        // simpler method assumes the correct allele and that peptide binds
        ++mTotalBinds;

        for(int pos = 0; pos < peptide.length(); ++pos)
        {
            char aminoAcid = peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            mBindCounts[aaIndex][pos] += 1;
        }
    }

    public double[][] getCounts(final String dataType)
    {
        if(dataType.equals(DATA_TYPE_BIND_COUNTS))
            return mBindCounts;
        else if(dataType.equals(DATA_TYPE_NOISE))
            return mNoiseCounts;
        if(dataType.equals(DATA_TYPE_LENGTH_WEIGHTED))
            return mWeightedCounts;
        else if(dataType.equals(DATA_TYPE_ALLELE_WEIGHTED))
            return mFinalWeightedCounts;
        else
            return null;
    }

    public static void writeCounts(final BufferedWriter writer, final BindCountData bindCounts, int maxPeptideLength, boolean writeNoise)
    {
        for(int i = 0; i < COUNT_DATA_TYPES.size(); ++i)
        {
            String dataType = COUNT_DATA_TYPES.get(i);

            if(!writeNoise && i == 1)
                continue;

            final double[][] counts = bindCounts.getCounts(dataType);

            writeCounts(writer, bindCounts, dataType, counts, maxPeptideLength);
        }
    }

    public static void writeCounts(
            final BufferedWriter writer, final BindCountData bindCounts, final String dataType, final double[][] counts, int maxPeptideLength)
    {
        try
        {
            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                writer.write(String.format("%s,%s,%d,%c", dataType, bindCounts.Allele, bindCounts.PeptideLength, aminoAcid));

                for(int pos = 0; pos < maxPeptideLength; ++pos)
                {
                    if(pos < bindCounts.PeptideLength)
                    {
                        writer.write(String.format(",%.1f", counts[aa][pos]));
                    }
                    else
                    {
                        writer.write(",0.0");
                    }
                }

                writer.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write bind counts data: {}", e.toString());
        }
    }

    public static BufferedWriter initFrequencyWriter(final String filename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,PeptideLength,AminoAcid,PeptidePos,BindCount");
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
        NE_LOGGER.debug("writing allele({}) frequency", Allele);

        try
        {
            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                for(int p = 0; p < PeptideLength; ++p)
                {
                    double totalBinds = mBindCounts[aa][p];

                    writer.write(String.format("%s,%d,%s,%d,%d", Allele, PeptideLength, aminoAcid, p, totalBinds));

                    writer.newLine();
                }
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write frequency data: {}", e.toString());
        }
    }

    public String toString()
    {
        return String.format("%s L%d: totalBinds(%d)", Allele, PeptideLength, mTotalBinds);
    }
}
