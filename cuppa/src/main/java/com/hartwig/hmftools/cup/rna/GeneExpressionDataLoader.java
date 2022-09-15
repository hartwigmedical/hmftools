package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
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

            CUP_LOGGER.info("loaded RNA expression data for {} samples from{}", sampleTpmIndex.size(), filename);

            return sampleGeneExpression;
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to load RNA expression ref data from {}: {}", filename, e.toString());
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
}
