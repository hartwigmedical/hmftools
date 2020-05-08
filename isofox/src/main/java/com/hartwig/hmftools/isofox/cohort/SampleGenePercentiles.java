package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.cohort.TransExpressionDistribution.DISTRIBUTION_SIZE;
import static com.hartwig.hmftools.isofox.cohort.TransExpressionDistribution.getTpmMedian;
import static com.hartwig.hmftools.isofox.cohort.TransExpressionDistribution.getTpmPercentile;
import static com.hartwig.hmftools.isofox.cohort.TransExpressionDistribution.loadCohortDistribution;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public class SampleGenePercentiles
{
    private final CohortConfig mConfig;

    private final Map<String,Map<String,double[]>> mCancerTypesGeneDistribution;

    private BufferedWriter mWriter;

    private static final String PAN_CANCER = "ALL";
    private static final int FLD_CANCER_TYPE = 0;
    private static final int FLD_FILENAME = 1;

    public SampleGenePercentiles(final CohortConfig config)
    {
        mConfig = config;
        mCancerTypesGeneDistribution = Maps.newHashMap();

        final Map<String,String> cancerGeneFilenames = loadCancerGeneDistributionFilenames(mConfig.CancerGeneFiles);

        for(Map.Entry<String,String> entry : cancerGeneFilenames.entrySet())
        {
            final String cancerType = entry.getKey();
            final String filename = entry.getValue();

            ISF_LOGGER.debug("cancerType({}) loading gene records", cancerType);

            final Map<String,double[]> geneDistributionMap = Maps.newHashMap();
            loadCohortDistribution(filename, geneDistributionMap, "gene", DISTRIBUTION_SIZE + 2, mConfig.RestrictedGeneIds);

            ISF_LOGGER.info("cancerType({}) loaded {} gene records", cancerType, geneDistributionMap.size());

            mCancerTypesGeneDistribution.put(cancerType, geneDistributionMap);
        }

        mWriter = null;
    }

    public void processSampleFiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, CohortAnalysisType.GENE_DISTRIBUTION, filenames))
            return;

        initialiseWriter();

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path genesFile = filenames.get(i);

            processSampleFile(sampleId, genesFile);
            ISF_LOGGER.debug("{}: sample({}) processed gene file", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples gene files", mConfig.SampleData.SampleIds.size());

        closeBufferedWriter(mWriter);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("sample_gene_perc_data.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,GeneId,GeneName,TPM,CohortMedianTPM,CancerMedianTPM,CohortPercentile,CancerPercentile");
            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write gene data file: {}", e.toString());
        }
    }

    private void writeSamplePercentileData(final String sampleId, final String geneId, final String geneName, double tpm)
    {
        try
        {
            final String cancerType = mConfig.SampleData.SampleCancerType.get(sampleId);

            try
            {
                double panCancerPerc = getTpmPercentile(mCancerTypesGeneDistribution.get(PAN_CANCER), geneId, tpm);
                double cancerPerc = getTpmPercentile(mCancerTypesGeneDistribution.get(cancerType), geneId, tpm);

                double panCancerMedian = getTpmMedian(mCancerTypesGeneDistribution.get(PAN_CANCER), geneId);
                double cancerMedian = getTpmMedian(mCancerTypesGeneDistribution.get(cancerType), geneId);

                mWriter.write(String.format("%s,%s,%s,%6.3e", sampleId, geneId, geneName, tpm));

                mWriter.write(String.format(",%6.3e,%6.3e,%.1f,%.1f", panCancerMedian, cancerMedian, panCancerPerc, cancerPerc));
                mWriter.newLine();
            }
            catch (NullPointerException e)
            {
                ISF_LOGGER.error("sampleId({}) geneId({}) cancerType({}) missing data", sampleId, geneId, cancerType, e.toString());
                return;
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write sample gene percentile data file: {}", e.toString());
        }
    }

    private void processSampleFile(final String sampleId, final Path filename)
    {
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

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(geneId))
                    continue;

                final String geneName = items[geneNameIndex];
                double tpm = Double.parseDouble(items[tpmIndex]);
                writeSamplePercentileData(sampleId, geneId, geneName, tpm);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
            return;
        }
    }

    private Map<String,String> loadCancerGeneDistributionFilenames(final String inputFile)
    {
        final Map<String,String> cancerGeneFileMap = Maps.newHashMap();

        if(!Files.exists(Paths.get(inputFile)))
        {
            ISF_LOGGER.error("invalid cancer-gene distribution file list");
            return cancerGeneFileMap;
        }

        try
        {
            final List<String> lines = Files.readAllLines(new File(inputFile).toPath());

            for(String data : lines)
            {
                final String[] items = data.split(",");
                final String cancerType = items[FLD_CANCER_TYPE];
                final String filename = items[FLD_FILENAME];

                cancerGeneFileMap.put(cancerType, filename);
            }

            ISF_LOGGER.info("loaded {} cancer-gene distribution filenames from file({})", cancerGeneFileMap.size(), inputFile);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load cancer-gene distribution file({}): {}", inputFile, e.toString());
        }

        return cancerGeneFileMap;
    }

}
