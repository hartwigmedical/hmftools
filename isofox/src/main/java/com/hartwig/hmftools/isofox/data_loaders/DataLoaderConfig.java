package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.IsofoxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.EXCLUDED_GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.IsofoxConfig.OUTPUT_ID;
import static com.hartwig.hmftools.isofox.IsofoxConfig.loadGeneIdsFile;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.getFileId;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class DataLoaderConfig
{
    public static final String ROOT_DATA_DIRECTORY = "root_data_dir";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    public static final String USE_SAMPLE_DIRS = "use_sample_dir";
    public static final String LOAD_TYPES = "load_types";
    public static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";

    public final String RootDataDir;
    public final String OutputDir;
    public final String OutputIdentifier;
    public final SampleDataCache SampleData;
    public final boolean UseSampleDirectories;

    public final int AltSJMinSampleThreshold;
    public final List<String> RestrictedGeneIds;
    public final List<String> ExcludedGeneIds;

    public final List<DataLoadType> LoadTypes;

    public DataLoaderConfig(final CommandLine cmd)
    {
        RootDataDir = cmd.getOptionValue(ROOT_DATA_DIRECTORY);
        UseSampleDirectories = cmd.hasOption(USE_SAMPLE_DIRS);

        String outputdir = cmd.getOptionValue(DATA_OUTPUT_DIR);
        if(!outputdir.endsWith(File.separator))
            outputdir += File.separator;
        OutputDir = outputdir;
        OutputIdentifier = cmd.getOptionValue(OUTPUT_ID);

        final String sampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE);

        SampleData = new SampleDataCache(sampleDataFile);

        if(!SampleData.isValid())
        {
            ISF_LOGGER.warn("invalid sample data file({})", sampleDataFile);
        }

        LoadTypes = Arrays.stream(cmd.getOptionValue(LOAD_TYPES).split(";"))
                .map(x -> DataLoadType.valueOf(x)).collect(Collectors.toList());

        AltSJMinSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_SAMPLES));

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            loadGeneIdsFile(inputFile, RestrictedGeneIds);
            ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }

        if(cmd.hasOption(EXCLUDED_GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(EXCLUDED_GENE_ID_FILE);
            loadGeneIdsFile(inputFile, ExcludedGeneIds);
            ISF_LOGGER.info("file({}) loaded {} excluded genes", inputFile, ExcludedGeneIds.size());
        }
    }

    public static boolean isValid(final CommandLine cmd)
    {
        return cmd.hasOption(ROOT_DATA_DIRECTORY) && cmd.hasOption(DATA_OUTPUT_DIR) && cmd.hasOption(SAMPLE_DATA_FILE)
                && cmd.hasOption(LOAD_TYPES);
    }

    public String formCohortFilename(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + "isfox_" + OutputIdentifier + "." + fileId;
        else
            return OutputDir + "isofox_" + fileId;
    }

    public static boolean formSampleFilenames(final DataLoaderConfig config, final DataLoadType dataType, final List<Path> filenames)
    {
        String rootDir = config.RootDataDir;

        if(!rootDir.endsWith(File.separator))
            rootDir += File.separator;

        for(final String sampleId : config.SampleData.SampleIds)
        {
            String filename = rootDir;

            if(config.UseSampleDirectories)
                filename += File.separator + sampleId + File.separator;

            filename += sampleId + ".isf.";
            filename += getFileId(dataType);

            final Path path = Paths.get(filename);

            if (!Files.exists(path))
            {
                ISF_LOGGER.error("sampleId({}) no file({}) found", sampleId, filename);
                filenames.clear();
                return false;
            }

            filenames.add(path);
        }

        return true;
    }


    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(ROOT_DATA_DIRECTORY, true, "Root data directory for input files or sample directories");
        options.addOption(SAMPLE_DATA_FILE, true, "File with list of samples and cancer types to load data for");
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(USE_SAMPLE_DIRS, false, "File with list of samples to load data for");
        options.addOption(LOAD_TYPES, true, "List of data types to load & process");
        options.addOption(ALT_SJ_MIN_SAMPLES, true, "Min numbner of samples to report an alt SJ");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(EXCLUDED_GENE_ID_FILE, true, "Optional CSV file of genes to ignore");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        return options;
    }
}
