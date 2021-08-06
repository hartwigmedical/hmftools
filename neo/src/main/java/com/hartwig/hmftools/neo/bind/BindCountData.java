package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class BindCountData
{
    public final String Allele;
    public final int PeptideLength;

    private final int[][] mObservations; // data points with a given amino acid and peptide
    private final double[][] mBindScoreTotals; // log-score total
    private final double[][] mBindCounts; // observations with affinity below configured binding threshold
    private final double[][] mNoiseCounts; // adjusted for noise

    private final double[][] mWeightedCounts; // weighted across the peptide lengths for this allele, padded to reference length
    private final double[][] mFinalWeightedCounts; // weighted across all other alleles, padded to reference length

    private final Map<String,double[]> mComboData;
    private int mTotal;
    private int mTotalBinds; // vs high affinity threshold

    public BindCountData(final String allele, final int peptideLength)
    {
        Allele = allele;
        PeptideLength = peptideLength;
        int aminoAcidCount = AMINO_ACID_COUNT;

        mObservations = new int[aminoAcidCount][PeptideLength];
        mBindScoreTotals = new double[aminoAcidCount][PeptideLength];
        mBindCounts = new double[aminoAcidCount][PeptideLength];
        mNoiseCounts = new double[aminoAcidCount][PeptideLength];
        mWeightedCounts = new double[aminoAcidCount][PeptideLength];
        mFinalWeightedCounts = new double[aminoAcidCount][PeptideLength];

        mComboData = Maps.newHashMap();
        mTotal = 0;
        mTotalBinds = 0;
    }

    public final double[][] getBindCounts() { return mBindCounts; }
    public final int[][] getObservations() { return mObservations; }
    public Map<String,double[]> getComboData() { return mComboData; }
    public final double[][] getNoiseCounts() { return mNoiseCounts; }
    public final double[][] getWeightedCounts() { return mWeightedCounts; }
    public final double[][] getFinalWeightedCounts() { return mFinalWeightedCounts; }
    public int totalBindCount() { return mTotalBinds; }

    public void logStats()
    {
        NE_LOGGER.info("allele({}) peptideLen({}) total({}) binds({}) pairs({})",
                Allele, PeptideLength, mTotal, mTotalBinds, mComboData.size());
    }

    public void processBindData(final BindData bindData, boolean calcPairs, final CalcConstants calcConstants)
    {
        if(bindData.Peptide.length() != PeptideLength || !bindData.Allele.equals(Allele))
            return;

        double levelScore = calcConstants.deriveLevelScore(bindData.Affinity);
        double bindPerc = calcConstants.deriveAffinityPercent(bindData.Affinity);

        ++mTotal;

        boolean actualBind = bindData.Affinity < calcConstants.BindingAffinityHigh;
        boolean predictedBind = bindData.predictedAffinity() < calcConstants.BindingAffinityHigh && actualBind;

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
            mBindCounts[aaIndex][pos] += bindPerc;

            if(calcPairs && bindData.isTraining())
            {
                if(pos < bindData.Peptide.length() - 1)
                {
                    for(int pos2 = pos + 1; pos2 < bindData.Peptide.length(); ++pos2)
                    {
                        char aminoAcid2 = bindData.Peptide.charAt(pos2);
                        if(aminoAcidIndex(aminoAcid2) == INVALID_AMINO_ACID)
                            continue;

                        ComboCorrelations.updatePairData(
                                mComboData, aminoAcid, pos, aminoAcid2, pos2, levelScore, actualBind, predictedBind);
                    }
                }
            }
        }
    }

    public void processBindingPeptide(final String peptide)
    {
        // simpler method assumes the correct allele and that peptide binds
        ++mTotal;
        ++mTotalBinds;

        for(int pos = 0; pos < peptide.length(); ++pos)
        {
            char aminoAcid = peptide.charAt(pos);
            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                continue;

            ++mObservations[aaIndex][pos];
            mBindScoreTotals[aaIndex][pos] += 1;
            mBindCounts[aaIndex][pos] += 1;
        }
    }

    private static final List<String> COUNT_DATA_TYPES = Lists.newArrayList(
            "BindCounts", "NoiseCounts", "PeptideLengthWeighted", "AlleleMotifWeighted");

    public static void writeCounts(final BufferedWriter writer, final BindCountData bindCounts, int maxPeptideLength, boolean writeNoise)
    {
        try
        {
            for(int i = 0; i < COUNT_DATA_TYPES.size(); ++i)
            {
                String dataType = COUNT_DATA_TYPES.get(i);

                if(!writeNoise && i == 1)
                    continue;

                final double[][] counts = (i == 0) ? bindCounts.getBindCounts() : (i == 1) ? bindCounts.getNoiseCounts() :
                        (i == 2) ? bindCounts.getWeightedCounts() : bindCounts.getFinalWeightedCounts();

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
        NE_LOGGER.debug("writing allele({}) frequency", Allele);

        try
        {
            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                for(int p = 0; p < PeptideLength; ++p)
                {
                    int freq = mObservations[aa][p];
                    double totalBinds = mBindCounts[aa][p];
                    double totalScore = mBindScoreTotals[aa][p];
                    double avgScore = freq > 0 ? totalScore / freq : 0;

                    writer.write(String.format("%s,%d,%s,%d,%d,%.4f,%.4f,%.2f",
                            Allele, PeptideLength, aminoAcid, p, freq, avgScore, totalScore, totalBinds));

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
