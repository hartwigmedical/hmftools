package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

public class GlobalWeights
{
    private final CalcConstants mConstants;

    private final Map<Integer,double[][]> mGlobalWeights; // sum of peptide-length weighted counts from all alleles
    private final Map<Integer,Double> mGlobalPeptideLengthTotals;

    public GlobalWeights(final CalcConstants calcConstants)
    {
        mConstants = calcConstants;

        mGlobalWeights = Maps.newHashMap();
        mGlobalPeptideLengthTotals = Maps.newHashMap();
    }

    public boolean enabled() { return mConstants.GlobalWeight > 0; }

    public void processBindCounts(final BindCountData bindCounts)
    {
        if(!enabled())
            return;

        final double[][] weightedCounts = bindCounts.getWeightedCounts();
        double[][] globalCounts = mGlobalWeights.get(bindCounts.PeptideLength);

        if(globalCounts == null)
        {
            globalCounts = new double[AMINO_ACID_COUNT][bindCounts.PeptideLength];
            mGlobalWeights.put(bindCounts.PeptideLength, globalCounts);
        }

        for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
        {
            for(int pos = 0; pos < bindCounts.PeptideLength; ++pos)
            {
                globalCounts[aa][pos] += weightedCounts[aa][pos];
            }
        }
    }

    public void setGlobalTotals()
    {
        if(!mGlobalPeptideLengthTotals.isEmpty() || !enabled())
            return;

        for(Map.Entry<Integer,double[][]> entry : mGlobalWeights.entrySet())
        {
            int peptideLength = entry.getKey();
            double[][] counts = entry.getValue();

            double total = Matrix.sumMatrix(counts);
            mGlobalPeptideLengthTotals.put(peptideLength, total);
        }
    }

    public void writeGlobalCounts(final BufferedWriter writer, int maxPeptideLength)
    {
        try
        {
            for(Map.Entry<Integer,double[][]> entry : mGlobalWeights.entrySet())
            {
                int peptideLength = entry.getKey();
                double[][] counts = entry.getValue();

                for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
                {
                    char aminoAcid = AMINO_ACIDS.get(aa);

                    writer.write(String.format("%s,%s,%d,%c", "Global", "GLOBAL", peptideLength, aminoAcid));

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
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write global counts data: {}", e.toString());
        }
    }

}
