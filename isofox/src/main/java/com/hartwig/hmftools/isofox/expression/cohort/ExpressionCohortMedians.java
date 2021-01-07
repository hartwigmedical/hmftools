package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.utils.MatrixUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_NAME;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class ExpressionCohortMedians
{
    // for re-writing expression matrix as cohort and cancer-type median values
    private final CohortConfig mConfig;
    private BufferedWriter mWriter;

    private Matrix mExpressionMatrix;
    private final Map<String,Integer> mSampleIndexMap;
    private final List<String> mGeneTransIds;
    private final Map<String,String> mIdNameMap;
    private boolean mTranscriptScope;

    public ExpressionCohortMedians(final CohortConfig config)
    {
        mConfig = config;

        mExpressionMatrix = null;
        mSampleIndexMap = Maps.newHashMap();
        mGeneTransIds = Lists.newArrayList();
        mIdNameMap = Maps.newHashMap();
        mTranscriptScope = false;

        mWriter = null;
    }

    public void produceCohortData()
    {
        if(!loadExpressionMatrix())
            return;

        initialiseWriter();

        for(int i = 0; i < mGeneTransIds.size(); ++i)
        {
            produceCohortData(i);
        }

        closeBufferedWriter(mWriter);
    }

    private void produceCohortData(int index)
    {
        final String geneTransId = mGeneTransIds.get(index);
        final String geneTransName = mIdNameMap.get(geneTransId);

        final double[] tpmValues = mExpressionMatrix.getRow(index);

        final List<Double> allValues = Lists.newArrayListWithExpectedSize(mExpressionMatrix.Cols);

        for(int s = 0; s < tpmValues.length; ++s)
        {
            addSortedTpm(allValues, tpmValues[s]);
        }

        double allMedian = calculatedMedian(allValues);

        try
        {
            if(mTranscriptScope)
                mWriter.write(String.format("%s,%s", geneTransName, geneTransId));
            else
                mWriter.write(String.format("%s,%s", geneTransId, geneTransName));

            mWriter.write(String.format(",%6.3e", allMedian));

            for(Map.Entry<String,List<String>> entry : mConfig.SampleData.CancerTypeSamples.entrySet())
            {
                final String cancerType = entry.getKey();
                final List<String> samples = entry.getValue();

                final List<Double> cancerValues = Lists.newArrayListWithExpectedSize(samples.size());

                for(final String sampleId : samples)
                {
                    Integer sampleIndex = mSampleIndexMap.get(sampleId);

                    if(sampleIndex == null)
                        continue;

                    addSortedTpm(cancerValues, tpmValues[sampleIndex]);
                }

                double cancerMedian = calculatedMedian(cancerValues);
                mWriter.write(String.format(",%6.3e", cancerMedian));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private double calculatedMedian(final List<Double> values)
    {
        if((values.size() % 2) == 0)
        {
            int medianIndexLow = values.size() / 2;
            int medianIndexHigh = medianIndexLow + 1;
            return (values.get(medianIndexLow) + values.get(medianIndexHigh)) * 0.5;
        }
        else
        {
            int medianIndex = values.size() / 2;
            return values.get(medianIndex);
        }
    }

    private void addSortedTpm(final List<Double> values, double tpm)
    {
        int index = 0;
        while(index < values.size())
        {
            if(tpm < values.get(index))
                break;

            ++index;
        }

        values.add(index, tpm);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mTranscriptScope ?
                    mConfig.formCohortFilename("transcript_medians.csv") : mConfig.formCohortFilename("gene_medians.csv");

            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("GeneId,GeneName");

            if(mTranscriptScope)
                mWriter.write(",TransName");

            for(String cancerType : mConfig.SampleData.CancerTypeSamples.keySet())
            {
                mWriter.write(String.format(",%s", cancerType));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
        }
    }

    private boolean loadExpressionMatrix()
    {
        // keep track of gene/transcript ids and names
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.Expression.GeneExpMatrixFile));

            String header = fileReader.readLine();

            final Map<String,Integer> fieldsMapIndex = createFieldsIndexMap(header, DELIMITER);

            int geneIdIndex = fieldsMapIndex.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMapIndex.get(FLD_GENE_NAME);
            int transNameIndex = fieldsMapIndex.containsKey(FLD_TRANS_NAME) ? fieldsMapIndex.get(FLD_TRANS_NAME) : -1;
            mTranscriptScope = transNameIndex >= 0;

            String line = fileReader.readLine();

            while(line != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String geneName = items[geneNameIndex];
                final String geneId = items[geneIdIndex];

                if(mTranscriptScope)
                {
                    final String transName = items[transNameIndex];
                    mGeneTransIds.add(transName);
                    mIdNameMap.put(transName, String.format("%s,%s", geneId, geneName));
                }
                else
                {
                    mGeneTransIds.add(geneId);
                    mIdNameMap.put(geneId, geneName);
                }

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            ISF_LOGGER.debug("failed to load RNA expression data from {}: {}", mConfig.Expression.GeneExpMatrixFile, e.toString());
            return false;
        }

        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_GENE_NAME);

        if(mTranscriptScope)
            ignoreFields.add(FLD_TRANS_NAME);

        mExpressionMatrix = loadMatrixDataFile(mConfig.Expression.GeneExpMatrixFile, mSampleIndexMap, ignoreFields);
        mExpressionMatrix.cacheTranspose();

        ISF_LOGGER.debug("loaded genes({}) and {} samples expression matrix data", mGeneTransIds.size(), mExpressionMatrix.Cols);
        return true;
    }

}
