package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.MatrixFile.loadMatrixDataFile;
import static com.hartwig.hmftools.common.stats.FdrCalcs.calculateFDRs;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
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
import com.hartwig.hmftools.common.stats.MannWhitneyUTest;
import com.hartwig.hmftools.common.stats.MwuResult;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.stats.PValueResult;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class ExpressionCohortCompare
{
    private final CohortConfig mConfig;

    private Matrix mGeneExpressionMatrix;
    private final Map<String,Integer> mSampleIndexMap;
    private final List<String> mGeneIds;
    private final Map<String,String> mGeneIdNameMap;

    private BufferedWriter mWriter;

    private static final String CANCER_TYPE_ALL = "ALL";
    private static final String CANCER_TYPE_OTHER = "Other";

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
        final List<String> ignoreFields = Lists.newArrayList(FLD_GENE_ID, FLD_GENE_NAME);
        mGeneExpressionMatrix = loadMatrixDataFile(mConfig.Expression.GeneExpMatrixFile, mSampleIndexMap, ignoreFields);
        mGeneExpressionMatrix.cacheTranspose();

        // keep track of gene ids and names
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.Expression.GeneExpMatrixFile));

            String header = fileReader.readLine();

            final Map<String,Integer> fieldsMapIndex = createFieldsIndexMap(header, DELIMITER);

            int geneIdIndex = fieldsMapIndex.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMapIndex.get(FLD_GENE_NAME);

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

        ISF_LOGGER.debug("loaded genes({}) and {} samples expression matrix data", mGeneIds.size(), mGeneExpressionMatrix.Cols);
    }

    public void runAnalysis()
    {
        if(mConfig.SampleData.CohortNames.size() != 2)
        {
            ISF_LOGGER.error("2 distinct cohorts required");
            return;
        }

        if(!mConfig.SampleData.CancerTypeSamples.isEmpty())
        {
            mConfig.SampleData.CancerTypeSamples.keySet().stream()
                    .filter(x -> !x.equals(CANCER_TYPE_OTHER))
                    .forEach(x -> analyseCohorts(x));
        }

        analyseCohorts(CANCER_TYPE_ALL);

        closeBufferedWriter(mWriter);
    }

    private void analyseCohorts(final String cancerType)
    {
        final String cohortA = mConfig.SampleData.CohortNames.get(0);
        final String cohortB = mConfig.SampleData.CohortNames.get(1);
        final List<Integer> cohortASampleIndices = Lists.newArrayList();
        final List<Integer> cohortBSampleIndices = Lists.newArrayList();

        for(final String sampleId : mConfig.SampleData.SampleIds)
        {
            if(!cancerType.equals(CANCER_TYPE_ALL) && !cancerType.equals(mConfig.SampleData.SampleCancerType.get(sampleId)))
                continue;

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

        int totalSampleCount = cohortASampleIndices.size() + cohortBSampleIndices.size();

        ISF_LOGGER.info("cancerType({}) totalSamples({}) cohortA({} samples={}) cohortB({} samples={})",
                cancerType, totalSampleCount, cohortA, cohortASampleIndices.size(), cohortB, cohortBSampleIndices.size());

        if(cohortASampleIndices.isEmpty() || cohortBSampleIndices.isEmpty())
            return;

        double[] cohortAValues = new double[cohortASampleIndices.size()];
        double[] cohortBValues = new double[cohortBSampleIndices.size()];

        // for each gene, extract the expression for each cohort
        final MannWhitneyUTest mww = new MannWhitneyUTest();
        final List<PValueResult> pValueResults = Lists.newArrayList();
        final Map<String,MwuResult> mwuResults = Maps.newHashMap();

        for(int geneIndex = 0; geneIndex < mGeneIds.size(); ++geneIndex)
        {
            final String geneId = mGeneIds.get(geneIndex);

            populateCohortValues(cohortAValues, geneIndex, cohortASampleIndices);
            populateCohortValues(cohortBValues, geneIndex, cohortBSampleIndices);

            final MwuResult result = mww.calculate(cohortAValues, cohortBValues);
            mwuResults.put(geneId, result);
            pValueResults.add(new PValueResult(geneId, result.PValue));
        }

        ISF_LOGGER.debug("calculating FDRs for {} results", pValueResults.size());
        calculateFDRs(pValueResults);

        for(int i = 0; i < pValueResults.size(); ++i)
        {
            final PValueResult pValue = pValueResults.get(i);
            final MwuResult result = mwuResults.get(pValue.Id);

            final String higherRankedCohort = result.rankAverage1() > result.rankAverage2(totalSampleCount) ? cohortA : cohortB;

            writeResults(cancerType, pValue.Id, pValue, higherRankedCohort);
        }
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

            mWriter.write("CancerType,GeneId,GeneName,PValue,QValue,TestRank,HigherCohort");
            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene expression comparison file: {}", e.toString());
        }
    }

    private void writeResults(final String cancerType, final String geneId, final PValueResult pValue, final String higherCohort)
    {
        try
        {
            final String geneName = mGeneIdNameMap.get(geneId);

            mWriter.write(String.format("%s,%s,%s,%g,%g,%d,%s",
                    cancerType, geneId, geneName, pValue.PValue, pValue.QValue, pValue.Rank, higherCohort));
            mWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write gene expression comparison file: {}", e.toString());
            return;
        }
    }

}
