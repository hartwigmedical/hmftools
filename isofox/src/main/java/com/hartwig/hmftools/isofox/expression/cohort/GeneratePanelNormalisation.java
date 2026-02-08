package com.hartwig.hmftools.isofox.expression.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_DIST_DESC;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_DIST_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.PANEL_TPM_NORMALISATION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.CohortGenePercentiles.PAN_CANCER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class GeneratePanelNormalisation
{
    private final CohortConfig mConfig;

    private final List<String> mGeneIds;
    private final List<String> mGeneNames;
    private final List<Double> mWgsCohortMedians;

    public static final String FLD_TPM_ADJUST_FACTOR = "AdjustmentFactor";

    private static final double MIN_WGS_MEDIAN_TPM = 0.01;

    public GeneratePanelNormalisation(final CohortConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;

        mGeneIds = Lists.newArrayList();
        mGeneNames = Lists.newArrayList();
        mWgsCohortMedians = Lists.newArrayList();

        loadGeneIdFile(configBuilder.getValue(GENE_ID_FILE));
        buildWgsCohortMedians(configBuilder.getValue(GENE_DIST_FILE));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(GENE_DIST_FILE, false, GENE_DIST_DESC);
        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);
    }

    private void loadGeneIdFile(final String geneIdFile)
    {
        try
        {
            List<String> fileContents = Files.readAllLines(new File(geneIdFile).toPath());
            String header = fileContents.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);
            fileContents.remove(0);

            for(String line : fileContents)
            {
                String[] values = line.split(CSV_DELIM, -1);
                mGeneIds.add(values[fieldsIndexMap.get(FLD_GENE_ID)]);
                mGeneNames.add(values[fieldsIndexMap.get(FLD_GENE_NAME)]);
            }

            ISF_LOGGER.info("loaded {} gene IDs from {}", mGeneIds.size(), geneIdFile);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to read file({}): {}", geneIdFile, e.toString());
        }
    }

    private void buildWgsCohortMedians(final String geneDistributionFile)
    {
        CohortGenePercentiles cohortPercentiles = new CohortGenePercentiles(geneDistributionFile);

        for(String geneId : mGeneIds)
        {
            double cohortMedian = cohortPercentiles.getTpmMedian(geneId, PAN_CANCER);
            mWgsCohortMedians.add(cohortMedian);
        }
    }

    public void processSamples()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, PANEL_TPM_NORMALISATION, filenames))
            return;

        ISF_LOGGER.info("processing {} samples gene expression files", mConfig.SampleData.SampleIds.size());

        Map<String,List<Double>> geneTpms = Maps.newHashMap();

        // load each sample's alt SJs and consolidate into a single list
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path sampleFile = filenames.get(i);

            processSampleFile(sampleFile, geneTpms);
            ISF_LOGGER.debug("{}: sample({}) processed file", i, sampleId);
        }

        ISF_LOGGER.info("writing normalisation file");

        writeNormalisation(geneTpms);
    }

    private void processSampleFile(final Path filename, final Map<String,List<Double>> geneTpms)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);
            String fileDelim = inferFileDelimiter(filename.toString());

            Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), fileDelim);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int tpmIndex = fieldsMap.get(FLD_ADJ_TPM);

            for(final String data : lines)
            {
                final String[] values = data.split(fileDelim, -1);

                final String geneId = values[geneIdIndex];

                if(!mGeneIds.contains(geneId))
                    continue;

                double tpm = Double.parseDouble(values[tpmIndex]);

                List<Double> tpms = geneTpms.get(geneId);

                if(tpms == null)
                {
                    tpms = Lists.newArrayList();
                    geneTpms.put(geneId, tpms);
                }

                tpms.add(tpm);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
        }
    }

    private void writeNormalisation(final Map<String,List<Double>> geneTpms)
    {
        try
        {
            final String fileType = "panel_gene_normalisation.csv";
            final String filename = mConfig.formCohortFilename(fileType);

            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("GeneId,GeneName,WgsMedian,PanelMedian," + FLD_TPM_ADJUST_FACTOR);
            writer.newLine();

            for(int i = 0; i < mGeneIds.size(); ++i)
            {
                String geneId = mGeneIds.get(i);
                String geneName = mGeneNames.get(i);

                List<Double> tpms = geneTpms.get(geneId);
                double panelMedian = Doubles.median(tpms);

                double wgsMedian = max(mWgsCohortMedians.get(i), MIN_WGS_MEDIAN_TPM);
                double adjustmentFactor = panelMedian / wgsMedian;


                writer.write(String.format("%s,%s,%4.3e,%4.3e,%.6f",
                        geneId, geneName, wgsMedian, panelMedian, adjustmentFactor));

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write panel gene normalisation file: {}", e.toString());
        }
    }
}
