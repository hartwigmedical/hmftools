package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_C_FREQ_ADJUST;
import static com.hartwig.hmftools.neo.bind.BindScoreMatrix.writeMatrixData;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class GlobalWeights
{
    private final CalcConstants mConstants;
    private final AminoAcidFrequency mAminoAcidFrequency;

    private final Map<Integer,double[][]> mWeightedCounts; // sum of peptide-length weighted counts from all alleles
    private final Map<Integer,Double> mPeptideLengthTotals;
    private final Map<Integer,BindScoreMatrix> mPeptideLengthMatrixMap;

    public static final String GLOBAL_COUNTS = "GLOBAL";

    public GlobalWeights(final CalcConstants calcConstants, final AminoAcidFrequency aminoAcidFrequency)
    {
        mConstants = calcConstants;
        mAminoAcidFrequency = aminoAcidFrequency;

        mWeightedCounts = Maps.newHashMap();
        mPeptideLengthTotals = Maps.newHashMap();
        mPeptideLengthMatrixMap = Maps.newHashMap();
    }

    public boolean enabled() { return mConstants.GlobalWeight > 0; }

    public final Map<Integer,BindScoreMatrix> getMatrixMap() { return mPeptideLengthMatrixMap; }

    public void processBindCounts(final BindCountData bindCounts)
    {
        if(!enabled())
            return;

        final double[][] weightedCounts = bindCounts.getWeightedCounts();
        double[][] globalCounts = mWeightedCounts.get(bindCounts.PeptideLength);

        if(globalCounts == null)
        {
            globalCounts = new double[AMINO_ACID_COUNT][bindCounts.PeptideLength];
            mWeightedCounts.put(bindCounts.PeptideLength, globalCounts);
        }

        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
        {
            for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
            {
                globalCounts[aa][pos] += weightedCounts[aa][pos];
            }
        }
    }

    public void createMatrixData()
    {
        for(Map.Entry<Integer,double[][]> entry : mWeightedCounts.entrySet())
        {
            int peptideLength = entry.getKey();
            double[][] weightedCounts = entry.getValue();

            BindScoreMatrix matrix = new BindScoreMatrix(GLOBAL_COUNTS, peptideLength);
            mPeptideLengthMatrixMap.put(peptideLength, matrix);

            final double[][] data = matrix.getBindScores();

            for(int pos = 0; pos < peptideLength; ++pos)
            {
                for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                {
                    double aaFrequency = mAminoAcidFrequency.getAminoAcidFrequency(aa);

                    // Peptide Weight = log(max(posWeight(aa,pos) * AaAdjust, 0.005 * posWeightTotal)/AaFreq,2)
                    double adjustedCount = weightedCounts[aa][pos];

                    if(aa == 'C')
                        adjustedCount *= AMINO_ACID_C_FREQ_ADJUST;

                    double posWeight = log(2, adjustedCount / aaFrequency);
                    data[aa][pos] = posWeight;
                }
            }
        }
    }

//    public double calcScore(final BindCountData bindCounts)
//    {
//
//    }

    public void buildMatrixData()
    {
        if(!mPeptideLengthTotals.isEmpty() || !enabled())
            return;

        for(Map.Entry<Integer,double[][]> entry : mWeightedCounts.entrySet())
        {
            int peptideLength = entry.getKey();
            double[][] counts = entry.getValue();

            double total = Matrix.sumMatrix(counts);
            mPeptideLengthTotals.put(peptideLength, total);
        }

        createMatrixData();
    }

    public void writeGlobalCounts(final BufferedWriter writer, int maxPeptideLength, boolean writePosWeightMatrix, boolean writeBindCounts)
    {
        if(writePosWeightMatrix)
            mPeptideLengthMatrixMap.values().forEach(x -> writeMatrixData(writer, x, maxPeptideLength));

        if(writeBindCounts)
        {
            try
            {
                for(Map.Entry<Integer, double[][]> entry : mWeightedCounts.entrySet())
                {
                    int peptideLength = entry.getKey();
                    double[][] counts = entry.getValue();

                    for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                    {
                        char aminoAcid = AMINO_ACIDS.get(aa);

                        writer.write(String.format("%s,%s,%d,%c", "Global", GLOBAL_COUNTS, peptideLength, aminoAcid));

                        for(int pos = 0; pos < maxPeptideLength; ++pos)
                        {
                            if(pos < peptideLength)
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
            catch(IOException e)
            {
                NE_LOGGER.error("failed to write global counts data: {}", e.toString());
            }
        }
    }

}
