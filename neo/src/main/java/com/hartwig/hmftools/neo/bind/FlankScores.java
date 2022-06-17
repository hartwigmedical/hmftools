package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_POS_WEIGHTS;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_DATA_TYPE;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_OBSERVED_AA_POS_FREQ;
import static com.hartwig.hmftools.neo.bind.FlankCounts.DOWN_1;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.FlankCounts.START_AMINO_ACID_ID;
import static com.hartwig.hmftools.neo.bind.FlankCounts.START_AMINO_ACID_FREQ;
import static com.hartwig.hmftools.neo.bind.FlankCounts.STOP_AMINO_ACID_FREQ;
import static com.hartwig.hmftools.neo.bind.FlankCounts.TOTAL_FLANK_AA_COUNT;
import static com.hartwig.hmftools.neo.bind.FlankCounts.UP_1;
import static com.hartwig.hmftools.neo.bind.FlankCounts.flankAminoAcidIndex;

import static org.apache.commons.math3.util.FastMath.log;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.neo.utils.AminoAcidFrequency;

public class FlankScores
{
    // by amino acid and position
    private final double[][] mPosWeights;
    private boolean mHasData;

    public FlankScores()
    {
        mPosWeights = new double[FLANK_AMINO_ACID_COUNT][TOTAL_FLANK_AA_COUNT];
        mHasData = false;
    }

    public final double[][] getPosWeights() { return mPosWeights; }

    public boolean hasData() { return mHasData; }

    public double calcScore(final String upFlank, final String downFlank)
    {
        if(!mHasData)
            return 0;

        if(upFlank.isEmpty() && downFlank.isEmpty())
            return 0;

        double score = 0;

        if(!upFlank.isEmpty())
        {
            int flankLength = upFlank.length();

            if(flankLength >= 3)
            {
                score += getPosWeight(0, upFlank.charAt(0));
                score += getPosWeight(1, upFlank.charAt(1));
                score += getPosWeight(2, upFlank.charAt(2));
            }
            else if(flankLength >= 2)
            {
                score += getPosWeight(1, upFlank.charAt(0));
                score += getPosWeight(2, upFlank.charAt(1));
            }
            else
            {
                score += getPosWeight(2, upFlank.charAt(0));
            }
        }
        else
        {
            score += getPosWeight(UP_1, START_AMINO_ACID_ID);
        }

        if(!downFlank.isEmpty())
        {
            int flankLength = downFlank.length();

            if(flankLength >= 3)
            {
                score += getPosWeight(3, downFlank.charAt(0));
                score += getPosWeight(4, downFlank.charAt(1));
                score += getPosWeight(5, downFlank.charAt(2));
            }
            else if(flankLength >= 2)
            {
                score += getPosWeight(3, downFlank.charAt(0));
                score += getPosWeight(4, downFlank.charAt(1));
            }
            else
            {
                score += getPosWeight(3, downFlank.charAt(0));
            }
        }
        else
        {
            score += getPosWeight(DOWN_1, STOP_AMINO_ACID);
        }

        return score;
    }

    private double getPosWeight(int position, char aminoAcid)
    {
        int aaIndex = flankAminoAcidIndex(aminoAcid);

        if(aaIndex == INVALID_AMINO_ACID)
            return 0;

        return mPosWeights[aaIndex][position];
    }

    public void createMatrix(final int[][] bindCounts)
    {
        AminoAcidFrequency aminoAcidFrequency = new AminoAcidFrequency();

        for(int pos = 0; pos < TOTAL_FLANK_AA_COUNT; ++pos)
        {
            double posTotalCount = 0;

            for(int aa = 0; aa < FLANK_AMINO_ACID_COUNT; ++aa)
            {
                posTotalCount += bindCounts[aa][pos];
            }

            posTotalCount = max(posTotalCount, 0.001); // for training data without observed counts in a peptide length

            for(int aa = 0; aa < FLANK_AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = FLANK_AMINO_ACIDS.get(aa);

                if(aminoAcid == START_AMINO_ACID_ID && pos != UP_1)
                {
                    mPosWeights[aa][pos] = 0;
                    continue;
                }
                else if(aminoAcid == STOP_AMINO_ACID && pos != DOWN_1)
                {
                    mPosWeights[aa][pos] = 0;
                    continue;
                }

                double aaFrequency;

                if(aminoAcid == START_AMINO_ACID_ID)
                    aaFrequency = START_AMINO_ACID_FREQ;
                else if(aminoAcid == STOP_AMINO_ACID)
                    aaFrequency = STOP_AMINO_ACID_FREQ;
                else
                    aaFrequency = aminoAcidFrequency.getAminoAcidFrequency(aa);

                // Peptide Weight = log(max(posWeight(aa,pos), 0.005 * posWeightTotal)/AaFreq,2)
                // Peptide Weight new = log(max(posWeight(aa,pos)/posWeightTotal, 0.005)/AaFreq,2) - now a % weight
                double bindPerc = bindCounts[aa][pos] / posTotalCount;

                // handle very low observation counts
                if(aminoAcid != START_AMINO_ACID_ID && aminoAcid != STOP_AMINO_ACID)
                {
                    bindPerc = max(bindPerc, MIN_OBSERVED_AA_POS_FREQ);
                }

                if(bindPerc > 0)
                {
                    double posWeight = log(2, bindPerc / aaFrequency);
                    mPosWeights[aa][pos] = posWeight;
                }
            }
        }

        mHasData = true;
    }

    public void writeData(final BufferedWriter writer)
    {
        try
        {
            // write pos weights
            for(int aa = 0; aa < FLANK_AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = FLANK_AMINO_ACIDS.get(aa);

                writer.write(String.format("%s,%c", DATA_TYPE_POS_WEIGHTS, aminoAcid));

                for(int pos = 0; pos < TOTAL_FLANK_AA_COUNT; ++pos)
                {
                    writer.write(String.format(",%.6f", mPosWeights[aa][pos]));
                }

                writer.newLine();
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write flank pos-weight data: {}", e.toString());
        }
    }

    public boolean loadPosWeights(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            int aaIndex = fieldsIndexMap.get(FLD_AMINO_ACID);
            int dataTypeIndex = fieldsIndexMap.get(FLD_DATA_TYPE);
            int peptideStartIndex = aaIndex + 1;

            for(String line : lines)
            {
                String[] items = line.split(DELIM, -1);

                if(!items[dataTypeIndex].equals(DATA_TYPE_POS_WEIGHTS))
                    continue;

                // DataType,AminoAcid,U3,U2,U1,D1,D2,D3

                char aminoAcid = items[aaIndex].charAt(0);
                int aminoAcidIndex = flankAminoAcidIndex(aminoAcid);

                int flankPos = 0;
                for(int i = peptideStartIndex; i < items.length; ++i, ++flankPos)
                {
                    double value = Double.parseDouble(items[i]);
                    mPosWeights[aminoAcidIndex][flankPos] = value;
                }
            }

            mHasData = true;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load flanking pos-weights data file: {}" ,e.toString());
            return false;
        }

        return true;
    }

}
