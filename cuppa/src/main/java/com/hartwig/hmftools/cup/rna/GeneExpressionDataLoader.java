package com.hartwig.hmftools.cup.rna;

import static java.lang.Math.log;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

public class GeneExpressionDataLoader
{
    public static final List<String> GENE_EXP_IGNORE_FIELDS = Lists.newArrayList("GeneId", "GeneName");

    public static Matrix loadGeneExpressionMatrix(
            final String filename, final Map<String,Integer> sampleTpmIndex, final List<String> sampleNames,
            final List<String> geneIds, final List<String> geneNames)
    {
        if(filename == null || filename.isEmpty())
            return null;

        try
        {
            // populate the gene info and sample names
            BufferedReader fileReader = createBufferedReader(filename);

            String header = fileReader.readLine();

            String[] columns = header.split(DATA_DELIM, -1);

            if(columns.length < 3 || !columns[0].equals(FLD_GENE_ID) || !columns[1].equals(FLD_GENE_NAME))
            {
                CUP_LOGGER.error("invalid gene expression file header");
                return null;
            }

            // assumes GeneId, GeneName, Sample1, Sample2 etc
            for(int i = 2; i < columns.length; ++i)
            {
                String sampleId = columns[i];
                sampleNames.add(sampleId);
                sampleTpmIndex.put(sampleId, i - 2);
            }

            String line = null;

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DATA_DELIM, -1);
                geneIds.add(values[0]);
                geneNames.add(values[1]);
            }

            // now load the actual matrix data
            fileReader = createBufferedReader(filename);
            fileReader.readLine(); // skip header

            Matrix sampleGeneExpression = new Matrix(sampleNames.size(), geneIds.size());
            final double[][] sampleData = sampleGeneExpression.getData();

            int geneIndex = 0;
            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DATA_DELIM, -1);

                int sampleIndex = 0;
                for(int i = 2; i < values.length; ++i)
                {
                    Double value = !values[i].isEmpty() ? Double.parseDouble(values[i]) : 0;
                    sampleData[sampleIndex][geneIndex] = value; // transposed
                    ++sampleIndex;
                }

                ++geneIndex;
            }

            return sampleGeneExpression;
        }
        catch (IOException e)
        {
            CUP_LOGGER.debug("failed to load RNA expression ref data from {}: {}", filename, e.toString());
        }

        return null;
    }


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
        Matrix sampleMatrix = loadMatrixDataFile(filename, sampleIndexMap, GENE_EXP_IGNORE_FIELDS, true);

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

            Matrix matrix = new Matrix(1, geneIdIndexMap.size());

            // GeneId,GeneName, etc AdjTPM

            int geneIdCol = fieldsIndexMap.get("GeneId");
            int adjTPM = fieldsIndexMap.get("AdjTPM");

            int unknownGeneCount = 0;

            for(String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                String geneId = items[geneIdCol];
                double adjTpm = Double.parseDouble(items[adjTPM]);
                Integer geneIdIndex = geneIdIndexMap.get(geneId);

                if(geneIdIndex == null)
                {
                    CUP_LOGGER.trace("unknown geneId({}) in sample file({})", geneId, filename);
                    ++unknownGeneCount;
                    continue;
                }

                double logTpm = log(adjTpm + 1);

                matrix.set(0, geneIdIndex, logTpm);
            }

            if(unknownGeneCount > 0)
            {
                CUP_LOGGER.warn("sample file({}) ignored {} unknown genes", filename, unknownGeneCount);
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
