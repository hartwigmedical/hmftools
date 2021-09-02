package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DATA_TYPE_POS_WEIGHTS;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_DATA_TYPE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE_LEN;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACIDS;
import static com.hartwig.hmftools.neo.bind.BindConstants.AMINO_ACID_COUNT;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class BindScoreMatrix
{
    public final String Allele;
    public final int PeptideLength;

    private final double[][] mPosWeights; // by amino acid and position

    private static final double INVALID_SCORE = -1000;

    public BindScoreMatrix(final String allele, final int peptideLength)
    {
        Allele = allele;
        PeptideLength = peptideLength;

        int aminoAcidCount = AMINO_ACIDS.size();

        mPosWeights = new double[aminoAcidCount][PeptideLength];
    }

    public final double[][] getBindScores() { return mPosWeights; }

    public double calcScore(final String peptide)
    {
        if(peptide.length() != PeptideLength)
            return INVALID_SCORE; // for now

        double score = 0;

        for(int i = 0; i < peptide.length(); ++i)
        {
            char aminoAcid = peptide.charAt(i);

            int aaIndex = aminoAcidIndex(aminoAcid);

            if(aaIndex == INVALID_AMINO_ACID)
                return INVALID_SCORE;

            double aaPosScore = mPosWeights[aaIndex][i];

            score += aaPosScore;
        }

        return score;
    }

    public static BufferedWriter initMatrixWriter(final String filename, int peptideLength, boolean writeCounts)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("DataType,Allele,PeptideLength,AminoAcid");

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

    public static void writeMatrixData(final BufferedWriter writer, final BindScoreMatrix matrix, int maxPeptideLength)
    {
        final double[][] data = matrix.getBindScores();

        try
        {
            for(int aa = 0; aa < AMINO_ACID_COUNT; ++aa)
            {
                char aminoAcid = AMINO_ACIDS.get(aa);

                writer.write(String.format("%s,%s,%d,%c", DATA_TYPE_POS_WEIGHTS, matrix.Allele, matrix.PeptideLength, aminoAcid));

                for(int pos = 0; pos < maxPeptideLength; ++pos)
                {
                    if(pos < matrix.PeptideLength)
                    {
                        writer.write(String.format(",%.6f", data[aa][pos]));
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
            NE_LOGGER.error("failed to write matrix data: {}", e.toString());
        }
    }

    public static List<BindScoreMatrix> loadFromCsv(final String filename)
    {
        List<BindScoreMatrix> matrixList = Lists.newArrayList();

        BindScoreMatrix currentMatrix = null;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int dataTypeIndex = fieldsIndexMap.get(FLD_DATA_TYPE);
            int peptideLenIndex = fieldsIndexMap.get(FLD_PEPTIDE_LEN);
            int aaIndex = fieldsIndexMap.get(FLD_AMINO_ACID);
            int peptideStartIndex = max(aaIndex, dataTypeIndex) + 1;

            Set<String> alleles = Sets.newHashSet();

            for(String line : lines)
            {
                String[] items = line.split(DELIM, -1);

                // Allele,PeptideLength,AminoAcid,P0,P1,P2,P3,P4,P5,P6,P7,P8
                //B4001,9,A,1.0286,-4.7395,0.4656,-0.1505,-0.3065,-0.0971,-0.3378,0.4915,-3.4338
                String dataType = items[dataTypeIndex];

                if(!dataType.equals(DATA_TYPE_POS_WEIGHTS))
                    continue;

                String allele = items[alleleIndex];
                alleles.add(allele);

                int peptideLength = Integer.parseInt(items[peptideLenIndex]);

                if(currentMatrix == null || !currentMatrix.Allele.equals(allele) || currentMatrix.PeptideLength != peptideLength)
                {
                    currentMatrix = new BindScoreMatrix(allele, peptideLength);
                    matrixList.add(currentMatrix);
                }

                char aminoAcid = items[aaIndex].charAt(0);
                int aminoAcidIndex = aminoAcidIndex(aminoAcid);

                int peptidePos = 0;
                for(int i = peptideStartIndex; i < items.length; ++i, ++peptidePos)
                {
                    double value = Double.parseDouble(items[i]);
                    currentMatrix.getBindScores()[aminoAcidIndex][peptidePos] = value;

                    if(peptidePos == currentMatrix.PeptideLength - 1)
                        break;
                }
            }

            NE_LOGGER.info("loading matrix data for {} alleles from {}", alleles.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to load pos-weights data file: {}" ,e.toString());
        }

        return matrixList;
    }
}
