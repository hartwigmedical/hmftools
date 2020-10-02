package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.FdrCalcs.calculateFDRs;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.stats.PValueResult;
import com.hartwig.hmftools.isofox.cohort.CohortAnalysisType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

public class ExpressionCohortCompare
{
    private final CohortConfig mConfig;

    private final Map<String,Map<String,List<Double>>> mCohortGeneExpDataMap;
    private final Map<String,String> mGeneIdNameMap;

    private BufferedWriter mWriter;

    public ExpressionCohortCompare(final CohortConfig config)
    {
        mConfig = config;

        mCohortGeneExpDataMap = Maps.newHashMap();
        mGeneIdNameMap = Maps.newHashMap();

        mWriter = null;
    }

    public void processSamples()
    {
        if(mConfig.SampleData.CohortNames.size() != 2)
        {
            ISF_LOGGER.error("2 distinct cohorts required");
            return;
        }

        for(String cohortName : mConfig.SampleData.CohortNames)
        {
            mCohortGeneExpDataMap.put(cohortName, Maps.newHashMap());
        }

        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, CohortAnalysisType.GENE_DISTRIBUTION, filenames))
            return;

        initialiseWriter();

        // load each sample's gene expression records and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final String cohortName = mConfig.SampleData.SampleCohort.get(sampleId);

            final Path genesFile = filenames.get(i);

            loadFile(genesFile, cohortName);
            ISF_LOGGER.debug("{}: sample({}) loaded genes", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples gene files", mConfig.SampleData.SampleIds.size());

        compareGeneDistributions();
        closeBufferedWriter(mWriter);
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
            ISF_LOGGER.error("failed to write gene expresion comparison file: {}", e.toString());
        }
    }

    private void loadFile(final Path filename, final String cohortName)
    {
        final Map<String,List<Double>> geneExpMap = mCohortGeneExpDataMap.get(cohortName);

        try
        {
            final List<String> lines = Files.readAllLines(filename);

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMap.get(FLD_GENE_NAME);
            int tpmIndex = fieldsMap.get(FLD_TPM);

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);

                final String geneId = items[geneIdIndex];

                final String geneName = mGeneIdNameMap.get(geneId);
                if(geneName == null)
                    mGeneIdNameMap.put(geneId, items[geneNameIndex]);

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                    continue;

                double tpm = Double.parseDouble(items[tpmIndex]);
                addGeneTpmData(geneExpMap, geneId, tpm);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private void addGeneTpmData(final Map<String,List<Double>> geneExpMap, final String geneId, double tpm)
    {
        if(tpm < mConfig.Expression.TpmLogThreshold)
            return;

        List<Double> tpmList = geneExpMap.get(geneId);
        if(tpmList == null)
        {
            geneExpMap.put(geneId, Lists.newArrayList(tpm));
            return;
        }

        int index = 0;
        while(index < tpmList.size())
        {
            if(tpm > tpmList.get(index))
                break;

            ++index;
        }

        tpmList.add(index, tpm);
    }

    private static final int MIN_SAMPLES = 5;

    private void compareGeneDistributions()
    {
        final String testCohort = mConfig.SampleData.CohortNames.get(1);
        final String controlCohort = mConfig.SampleData.CohortNames.get(0);

        final Map<String,List<Double>> testGeneExpData = mCohortGeneExpDataMap.get(testCohort);
        final Map<String,List<Double>> controlGeneExpData = mCohortGeneExpDataMap.get(controlCohort);

        MannWhitneyUTest mww = new MannWhitneyUTest();

        final List<PValueResult> pValueResults = Lists.newArrayList();

        for(Map.Entry<String,List<Double>> entry : testGeneExpData.entrySet())
        {
            final String geneId = entry.getKey();
            final List<Double> testGeneExp = entry.getValue();
            final List<Double> controlGeneExp = controlGeneExpData.get(geneId);

            if(controlGeneExp == null || testGeneExp.size() < MIN_SAMPLES || controlGeneExp.size() < MIN_SAMPLES)
                continue;

            double[] testExpValues = convertList(testGeneExp);
            double[] controlExpValues = convertList(controlGeneExp);

            double pValue = mww.mannWhitneyUTest(testExpValues, controlExpValues);

            pValueResults.add(new PValueResult(geneId, pValue));
        }

        ISF_LOGGER.debug("calculating FDRs for {} results", pValueResults.size());
        calculateFDRs(pValueResults);

        for(final PValueResult pValue : pValueResults)
        {
            writeResults(pValue.Id, pValue);
        }
    }

    private void writeResults(final String geneId, final PValueResult pValue)
    {
        try
        {
            final String geneName = mGeneIdNameMap.get(geneId);

            mWriter.write(String.format("%s,%s,%g,%g,%d", geneId, geneName, pValue.PValue, pValue.QValue, pValue.Rank));
            mWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write gene expression comparison file: {}", e.toString());
            return;
        }
    }

}
