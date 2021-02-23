package com.hartwig.hmftools.isofox.loader;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.IsofoxConfig.loadGeneIdsFile;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class DataLoaderConfig
{
    public static final String SAMPLE_DATA_DIRECTORY = "sample_data_dir";
    public static final String SAMPLE = "sample";
    public static final String SAMPLE_ID_FILE = "sample_data_file";
    public static final String CANCER_TYPES_FILE = "cancer_types_file";
    public static final String GENE_DIST_FILE = "gene_distribution_file";

    public final String SampleDataDir;
    public final String GeneDistributionFile;
    public final List<String> RestrictedGeneIds;
    public final List<String> SampleIds;
    public final List<String> PrimaryCancerTypes;
    public final Map<String,String> SampleCancerTypes;

    public final DatabaseAccess DbAccess;

    public DataLoaderConfig(final CommandLine cmd)
    {
        SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIRECTORY));
        GeneDistributionFile = cmd.getOptionValue(GENE_DIST_FILE);

        SampleIds = Lists.newArrayList();
        PrimaryCancerTypes = Lists.newArrayList();
        SampleCancerTypes = Maps.newHashMap();

        if(cmd.hasOption(SAMPLE))
        {
            addSampleInfo(cmd.getOptionValue(SAMPLE));
        }
        else if(cmd.hasOption(SAMPLE_ID_FILE))
        {
            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(cmd.getOptionValue(SAMPLE_ID_FILE)));
                lines.stream().filter(x -> !x.equals("SampleId")).forEach(x -> addSampleInfo(x));
            }
            catch(IOException e)
            {
                ISF_LOGGER.warn("invalid sampleId file: {}", e.toString());
            }
        }

        if(cmd.hasOption(CANCER_TYPES_FILE))
        {
            try
            {
                final List<String> lines = Files.readAllLines(Paths.get(cmd.getOptionValue(CANCER_TYPES_FILE)));
                lines.stream().filter(x -> !x.equals("CancerType")).forEach(x -> PrimaryCancerTypes.add(x));
            }
            catch(IOException e)
            {
                ISF_LOGGER.warn("invalid sampleId file: {}", e.toString());
            }
        }

        RestrictedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            loadGeneIdsFile(inputFile, RestrictedGeneIds);
            ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }

        DbAccess = createDatabaseAccess(cmd);
    }

    private void addSampleInfo(final String sampleInfo)
    {
        if(!sampleInfo.contains(DELIMITER))
        {
            SampleIds.add(sampleInfo);
            return;
        }

        String[] items = sampleInfo.split(DELIMITER, -1);
        String sampleId = items[0];
        String cancerType = items[1];
        SampleIds.add(sampleId);
        SampleCancerTypes.put(sampleId, cancerType);
    }

    public boolean processGeneId(final String geneId)
    {
        return RestrictedGeneIds.isEmpty() || RestrictedGeneIds.contains(geneId);
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Sample tumor ID");
        options.addOption(SAMPLE_ID_FILE, true, "File with list of samples and cancer types to load data for");
        options.addOption(SAMPLE_DATA_DIRECTORY, true, "Root data directory for input files or sample directories");
        options.addOption(CANCER_TYPES_FILE, true, "Primary cancer types (otherwise will use 'Other' for sample");
        options.addOption(GENE_DIST_FILE, true, "Gene distribution for medians and percentile data");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
