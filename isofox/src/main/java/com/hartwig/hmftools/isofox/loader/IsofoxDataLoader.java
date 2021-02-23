package com.hartwig.hmftools.isofox.loader;

import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.isValid;
import static com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles.PAN_CANCER;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.GENE_RESULTS_FILE;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles;
import com.hartwig.hmftools.patientdb.dao.IsofoxDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Result;

public class IsofoxDataLoader
{
    private final DataLoaderConfig mConfig;
    private final SampleGenePercentiles mGeneDistribution;

    public IsofoxDataLoader(final CommandLine cmd)
    {
        mConfig = new DataLoaderConfig(cmd);
        mGeneDistribution = new SampleGenePercentiles(mConfig.GeneDistributionFile);
    }

    public boolean load()
    {
        if(mConfig.SampleIds.isEmpty())
        {
            ISF_LOGGER.error("no sample IDs configured");
            return false;
        }

        if(mConfig.DbAccess == null)
        {
            ISF_LOGGER.error("invalid DB connection");
            return false;
        }

        final IsofoxDAO rnaDAO = new IsofoxDAO(mConfig.DbAccess.context());

        for(String sampleId : mConfig.SampleIds)
        {
            final String sampleDataDir = mConfig.SampleDataDir.contains("*") ?
                    mConfig.SampleDataDir.replaceAll("\\*", sampleId) : mConfig.SampleDataDir;

            final String cancerType = getSampleCancerType(sampleId);

            ISF_LOGGER.debug("sample({}) cancerType({}) loading Isofox RNA data from {}", sampleId, cancerType, sampleDataDir);

            loadGeneExpression(sampleId, cancerType, sampleDataDir, rnaDAO);
        }

        return true;
    }

    private String getSampleCancerType(final String sampleId)
    {
        if(mConfig.SampleCancerTypes.containsKey(sampleId))
            return mConfig.SampleCancerTypes.get(sampleId);

        try
        {
            final String queryStr = String.format("select primaryTumorLocation from clinical where sampleId = '%s';", sampleId);

            final Result<Record> result = mConfig.DbAccess.context().fetch(queryStr);

            for(Record record : result)
            {
                return record.getValue("primaryTumorLocation").toString();
            }
        }
        catch(Exception e)
        {
            ISF_LOGGER.error("failed to retrieve clinical data for sample({}): {}", e.toString());
        }

        return CANCER_TYPE_OTHER;
    }

    private void loadGeneExpression(final String sampleId, final String cancerType, final String sampleDataDir, final IsofoxDAO rnaDAO)
    {
        final List<GeneExpression> geneExpressions = Lists.newArrayList();
        final String filename = sampleDataDir + sampleId + ISF_FILE_ID + GENE_RESULTS_FILE;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String geneId = items[fieldsMap.get(FLD_GENE_ID)];

                if(!mConfig.processGeneId(geneId))
                    continue;

                double tpm = Double.parseDouble(items[fieldsMap.get(FLD_TPM)]);

                double medianCancer = mGeneDistribution.getTpmMedian(geneId, cancerType);
                double percentileCancer = mGeneDistribution.getTpmPercentile(geneId, cancerType, tpm);
                double medianCohort = mGeneDistribution.getTpmMedian(geneId, PAN_CANCER);
                double percentileCohort = mGeneDistribution.getTpmPercentile(geneId, PAN_CANCER, tpm);

                geneExpressions.add(ImmutableGeneExpression.builder()
                        .geneName(items[fieldsMap.get(FLD_GENE_NAME)])
                        .tpm(tpm)
                        .splicedFragments(Integer.parseInt(items[fieldsMap.get(FLD_SPLICED_FRAGS)]))
                        .unsplicedFragments(Integer.parseInt(items[fieldsMap.get(FLD_UNSPLICED_FRAGS)]))
                        .medianTpmCancer(medianCancer)
                        .percentileCancer(percentileCancer)
                        .medianTpmCohort(medianCohort)
                        .percentileCohort(percentileCohort)
                        .build());
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename.toString(), e.toString());
            return;
        }

        ISF_LOGGER.debug("sample({}) writing {} gene expression records to DB", sampleId, geneExpressions.size());
        rnaDAO.writeGeneExpressions(sampleId, geneExpressions);
    }

    private void loadFusions(final String sampleId)
    {

    }


    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = DataLoaderConfig.createCmdLineOptions();
        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        IsofoxDataLoader dataLoader = new IsofoxDataLoader(cmd);

        if(!dataLoader.load())
        {
            ISF_LOGGER.info("Isofox data loading failed");
            return;
        }

        ISF_LOGGER.info("Isofox data loading complete");
    }

}
