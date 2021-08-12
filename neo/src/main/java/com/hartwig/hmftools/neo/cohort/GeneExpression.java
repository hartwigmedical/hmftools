package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.DELIMITER;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;

public class GeneExpression
{
    private final Matrix mGeneExpressionMatrix;
    private final Map<String,Integer> mSampleGeneExpIndexMap;
    private final Map<String,Integer> mGeneIdIndices;

    public GeneExpression(final String geneExpressionFile)
    {
        if(geneExpressionFile != null)
        {
            mSampleGeneExpIndexMap = Maps.newHashMap();
            mGeneIdIndices = Maps.newHashMap();
            mGeneExpressionMatrix = loadCohortGeneExpression(geneExpressionFile);
        }
        else
        {
            mSampleGeneExpIndexMap = null;
            mGeneIdIndices = null;
            mGeneExpressionMatrix = null;
        }
    }

    public double getExpression(final String geneId, final String sampleId)
    {
        Integer geneIndex = mGeneIdIndices.get(geneId);
        Integer sampleIndex = mSampleGeneExpIndexMap.get(sampleId);

        if(geneIndex == null || sampleIndex == null)
            return 0;

        return mGeneExpressionMatrix.get(geneIndex, sampleIndex);
    }

    private Matrix loadCohortGeneExpression(final String filename)
    {
        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_GENE_NAME);
        Matrix geneExpMatrix = loadMatrixDataFile(filename, mSampleGeneExpIndexMap, ignoreFields);

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            String header = fileData.get(0);
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
            fileData.remove(0);

            // GeneId,GeneName,CancerTypes..

            int geneIdCol = fieldsIndexMap.get(FLD_GENE_ID);

            for(int i = 0; i < fileData.size(); ++i)
            {
                final String[] items = fileData.get(i).split(DELIMITER, -1);
                final String geneId = items[geneIdCol];
                mGeneIdIndices.put(geneId, i);
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to read RNA gene Ids from file({}): {}", filename, e.toString());
            return null;
        }

        return geneExpMatrix;
    }

}
