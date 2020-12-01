package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.sigs.SigUtils.loadMatrixDataFile;
import static com.hartwig.hmftools.common.stats.FdrCalcs.calculateFDRs;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.SigMatrix;
import com.hartwig.hmftools.common.stats.PValueResult;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

public class ExpressionCohortCompare
{
    private final CohortConfig mConfig;

    private SigMatrix mGeneExpressionMatrix;
    private final Map<String,Integer> mSampleIndexMap;
    private final List<String> mGeneIds;
    private final Map<String,String> mGeneIdNameMap;

    private BufferedWriter mWriter;

    public ExpressionCohortCompare(final CohortConfig config)
    {
        mConfig = config;

        mSampleIndexMap = Maps.newHashMap();
        mGeneExpressionMatrix = null;

        mGeneIds = Lists.newArrayList();
        mGeneIdNameMap = Maps.newHashMap();

        loadGeneExpression();

        mWriter = null;
        initialiseWriter();
    }

    private void loadGeneExpression()
    {
        final List<String> ignoreFields = Lists.newArrayList("GeneId", "GeneName");
        mGeneExpressionMatrix = loadMatrixDataFile(mConfig.Expression.GeneExpMatrixFile, mSampleIndexMap, ignoreFields);
        mGeneExpressionMatrix.cacheTranspose();

        // keep track of gene ids and names
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.Expression.GeneExpMatrixFile));

            String header = fileReader.readLine();

            final Map<String,Integer> fieldsMapIndex = createFieldsIndexMap(header, DELIMITER);

            int geneIdIndex = fieldsMapIndex.get("GeneId");
            int geneNameIndex = fieldsMapIndex.get("GeneName");

            String line = fileReader.readLine();

            while(line != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                mGeneIds.add(items[geneIdIndex]);
                mGeneIdNameMap.put(items[geneIdIndex], items[geneNameIndex]);

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            ISF_LOGGER.debug("failed to load RNA expression ref data from {}: {}", mConfig.Expression.GeneExpMatrixFile, e.toString());
        }

        ISF_LOGGER.debug("loaded genes({}) and {} samples({}) expression matrix data", mGeneIds.size(), mGeneExpressionMatrix.Cols);
    }

    public void runAnalysis()
    {
        if(mConfig.SampleData.CohortNames.size() != 2)
        {
            ISF_LOGGER.error("2 distinct cohorts required");
            return;
        }

        final String cohortA = mConfig.SampleData.CohortNames.get(0);
        final String cohortB = mConfig.SampleData.CohortNames.get(1);
        final List<Integer> cohortASampleIndices = Lists.newArrayList();
        final List<Integer> cohortBSampleIndices = Lists.newArrayList();

        for(final String sampleId : mConfig.SampleData.SampleIds)
        {
            final String cohortName = mConfig.SampleData.SampleCohort.get(sampleId);
            Integer matrixIndex = mSampleIndexMap.get(sampleId);

            if(matrixIndex == null)
            {
                ISF_LOGGER.warn("missing sample({}) matrix data", sampleId);
                continue;
            }

            if(cohortName.equals(cohortA))
                cohortASampleIndices.add(matrixIndex);
            else
                cohortBSampleIndices.add(matrixIndex);
        }

        ISF_LOGGER.info("cohortA({} samples={}) cohortB({} samples={})",
                cohortA, cohortASampleIndices.size(), cohortB, cohortBSampleIndices.size());

        double[] cohortAValues = new double[cohortASampleIndices.size()];
        double[] cohortBValues = new double[cohortBSampleIndices.size()];

        // for each gene, extract the expression for each cohort
        final MannWhitneyUTest mww = new MannWhitneyUTest();
        final List<PValueResult> pValueResults = Lists.newArrayList();

        for(int geneIndex = 0; geneIndex < mGeneIds.size(); ++geneIndex)
        {
            final String geneId = mGeneIds.get(geneIndex);

            populateCohortValues(cohortAValues, geneIndex, cohortASampleIndices);
            populateCohortValues(cohortBValues, geneIndex, cohortBSampleIndices);

            double pValue = mww.mannWhitneyUTest(cohortAValues, cohortBValues);
            pValueResults.add(new PValueResult(geneId, pValue));
        }

        ISF_LOGGER.debug("calculating FDRs for {} results", pValueResults.size());
        calculateFDRs(pValueResults);

        for(final PValueResult pValue : pValueResults)
        {
            writeResults(pValue.Id, pValue);
        }

        closeBufferedWriter(mWriter);
    }

    private void populateCohortValues(final double[] cohortValues, int geneIndex, final List<Integer> cohortSampleIndices)
    {
        final double[][] matrixData = mGeneExpressionMatrix.getData();
        int index = 0;
        for(Integer sampleIndex : cohortSampleIndices)
        {
            cohortValues[index++] = matrixData[geneIndex][sampleIndex];
        }
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("gene_expression_compare.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("GeneId,GeneName,PValue,QValue,TestRank");
            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene expression comparison file: {}", e.toString());
        }
    }

    private void writeResults(final String geneId, final PValueResult pValue)
    {
        try
        {
            final String geneName = mGeneIdNameMap.get(geneId);

            mWriter.write(String.format("%s,%s,%g,%g,%d",
                    geneId, geneName, pValue.PValue, pValue.QValue, pValue.Rank));
            mWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write gene expression comparison file: {}", e.toString());
            return;
        }
    }

}
