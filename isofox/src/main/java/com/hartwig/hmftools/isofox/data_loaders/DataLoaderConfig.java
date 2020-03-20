package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.IsofoxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.IsofoxConfig.OUTPUT_ID;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class DataLoaderConfig
{
    public static final String ROOT_DATA_DIRECTORY = "root_data_dir";
    public static final String SAMPLE_IDS_FILE = "sample_ids_file";
    public static final String USE_SAMPLE_DIRS = "use_sample_dir";
    public static final String LOAD_TYPES = "load_types";

    public final String RootDataDir;
    public final String OutputDir;
    public final String OutputIdentifier;
    public final List<String> SampleIds;
    public final boolean UseSampleDirectories;

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

        SampleIds = Lists.newArrayList();
        final String sampleIdsFile = cmd.getOptionValue(SAMPLE_IDS_FILE);

        try
        {
            SampleIds.addAll(Files.readAllLines(new File(sampleIdsFile).toPath()));
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load gene ID file({}): {}", sampleIdsFile, e.toString());
        }

        LoadTypes = Arrays.stream(cmd.getOptionValue(LOAD_TYPES).split(";"))
                .map(x -> DataLoadType.valueOf(x)).collect(Collectors.toList());
    }

    public String formCohortFilename(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + "isfox_" + OutputIdentifier + "." + fileId;
        else
            return OutputDir + "isofox_" + fileId;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(ROOT_DATA_DIRECTORY, true, "Root data directory for input files or sample directories");
        options.addOption(SAMPLE_IDS_FILE, true, "File with list of samples to load data for");
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(USE_SAMPLE_DIRS, false, "File with list of samples to load data for");
        options.addOption(LOAD_TYPES, true, "List of data types to load & process");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        return options;
    }
}
