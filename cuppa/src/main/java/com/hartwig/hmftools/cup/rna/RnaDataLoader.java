package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.log;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.utils.Matrix;

public class RnaDataLoader
{
    public static final List<String> GENE_EXP_IGNORE_FIELDS = Lists.newArrayList("GeneId", "GeneName");

    public static boolean loadGeneIdIndices(final String filename, final Map<String,Integer> geneIdIndices)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            String header = fileData.get(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            fileData.remove(0);

            // GeneId,GeneName,CancerTypes..

            int geneIdCol = fieldsIndexMap.get("GeneId");

            for(int i = 0; i < fileData.size(); ++i)
            {
                final String[] items = fileData.get(i).split(DATA_DELIM, -1);
                final String geneId = items[geneIdCol];
                geneIdIndices.put(geneId, i);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA gene Ids from file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static Matrix loadSampleGeneExpressionMatrix(
            final String filename, final Map<String,Integer> refGeneIdIndexMap, final Map<String,Integer> sampleIndexMap)
    {
        Matrix sampleMatrix = loadMatrixDataFile(filename, sampleIndexMap, GENE_EXP_IGNORE_FIELDS);

        // ensure genes are ordered as per the reference data
        final Map<String,Integer> sampleGeneIdIndices = Maps.newHashMap();
        if(!loadGeneIdIndices(filename, sampleGeneIdIndices))
            return null;

        if(sampleGeneIdIndices.size() != refGeneIdIndexMap.size())
        {
            CUP_LOGGER.error("sample gene expression matrix size({}) differs from ref size({})",
                    sampleGeneIdIndices.size(), refGeneIdIndexMap.size());
            return null;
        }

        if(sampleGeneIdIndices.entrySet().stream().anyMatch(x -> !refGeneIdIndexMap.containsKey(x.getKey())))
        {
            CUP_LOGGER.error("sample gene expression matrix missing ref gene entries");
            return null;
        }

        if(sampleGeneIdIndices.entrySet().stream().anyMatch(x -> refGeneIdIndexMap.get(x.getKey()) != x.getValue()))
        {
            // sample matrix is ordered differently so build it to match
            Matrix newMatrix = new Matrix(sampleMatrix.Rows, sampleMatrix.Cols);

            for(Map.Entry<String,Integer> refEntry : refGeneIdIndexMap.entrySet())
            {
                String refGeneId = refEntry.getKey();
                int refGeneIndex = refEntry.getValue();
                int sampleGeneIndex = sampleGeneIdIndices.get(refGeneId);

                newMatrix.setRow(refGeneIndex, sampleMatrix.getRow(sampleGeneIndex));
            }

            return newMatrix;
        }

        return sampleMatrix;
    }

    public static Matrix loadSampleGeneExpressionFile(final String filename, final Map<String,Integer> geneIdIndexMap)
    {
        if(!Files.exists(Paths.get(filename)))
            return null;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            String header = fileData.get(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            fileData.remove(0);

            Matrix matrix = new Matrix(geneIdIndexMap.size(), 1);

            // GeneId,GeneName, etc AdjTPM

            int geneIdCol = fieldsIndexMap.get("GeneId");
            int adjTPM = fieldsIndexMap.get("AdjTPM");

            for(String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                String geneId = items[geneIdCol];
                double adjTpm = Double.parseDouble(items[adjTPM]);
                Integer geneIdIndex = geneIdIndexMap.get(geneId);

                if(geneIdIndex == null)
                {
                    CUP_LOGGER.warn("unknown geneId({}) in sample file({})", geneId, filename);
                    continue;
                }

                double logTpm = log(adjTpm + 1);

                matrix.set(geneIdIndex, 0, logTpm);
            }

            return matrix;
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA sample gene data file({}): {}", filename, e.toString());
            return null;
        }
    }

}
