package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.isofox.IsofoxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.EXCLUDED_GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.IsofoxConfig.OUTPUT_ID;
import static com.hartwig.hmftools.isofox.IsofoxConfig.loadGeneIdsFile;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.getFileId;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ISOFOX_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.fusion.cohort.FusionCohortConfig;
import com.hartwig.hmftools.isofox.novel.cohort.AltSpliceJunctionCohort;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class CohortConfig
{
    public static final String ROOT_DATA_DIRECTORY = "root_data_dir";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    public static final String USE_SAMPLE_DIRS = "use_sample_dir";
    public static final String ALL_AVAILABLE_FILES = "all_available_files";
    public static final String LOAD_TYPES = "load_types";
    public static final String FAIL_MISSING = "fail_on_missing_file";

    public static final String TPM_LOG_THRESHOLD = "tpm_log_threshold";
    public static final String TPM_ROUNDING = "tpm_rounding";

    public static final String WRITE_SAMPLE_GENE_DISTRIBUTION_DATA = "write_sample_gene_dist";


    public static final String CANCER_GENE_FILES = "cancer_gene_files";

    public static final String COHORT_TRANS_FILE = "cohort_trans_file";
    public static final String CANCER_TRANS_FILE = "cancer_trans_file";

    public static final String SAMPLE_MUT_FILE = "sample_mut_file";
    private static final String THREADS = "threads";

    public final String RootDataDir;
    public final String OutputDir;
    public final String OutputIdentifier;
    public final SampleDataCache SampleData;
    public final boolean UseSampleDirectories;
    public final boolean AllAvailableFiles;
    public final List<String> RestrictedGeneIds;
    public final List<String> ExcludedGeneIds;
    public final boolean FailOnMissingSample;
    public final boolean ConvertUnmatchedCancerToOther;

    public final List<CohortAnalysisType> LoadTypes;

    // routine specifics
    public final String EnsemblDataCache;

    public final boolean WriteSampleGeneDistributionData;

    public final String CohortTransFile;
    public final String CancerTransFile;
    public final String SampleMutationsFile;

    public final String CancerGeneFiles;
    public final double TpmLogThreshold;
    public final double TpmRounding;

    public final FusionCohortConfig Fusions;

    public final int Threads;

    public CohortConfig(final CommandLine cmd)
    {
        String rootDir = cmd.getOptionValue(ROOT_DATA_DIRECTORY);

        if(!rootDir.endsWith(File.separator))
            rootDir += File.separator;

        RootDataDir = rootDir;

        UseSampleDirectories = cmd.hasOption(USE_SAMPLE_DIRS);
        AllAvailableFiles = !UseSampleDirectories && cmd.hasOption(ALL_AVAILABLE_FILES);
        FailOnMissingSample = cmd.hasOption(FAIL_MISSING);

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

        LoadTypes = Arrays.stream(cmd.getOptionValue(LOAD_TYPES).split(ITEM_DELIM))
                .map(x -> CohortAnalysisType.valueOf(x)).collect(Collectors.toList());

        CancerGeneFiles = cmd.getOptionValue(CANCER_GENE_FILES);
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

        EnsemblDataCache = cmd.getOptionValue(GENE_TRANSCRIPTS_DIR);

        ConvertUnmatchedCancerToOther = true;

        TpmLogThreshold = Double.parseDouble(cmd.getOptionValue(TPM_LOG_THRESHOLD, "0"));
        TpmRounding = Double.parseDouble(cmd.getOptionValue(TPM_ROUNDING, "2"));

        WriteSampleGeneDistributionData = cmd.hasOption(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA);
        CohortTransFile = cmd.getOptionValue(COHORT_TRANS_FILE);
        CancerTransFile = cmd.getOptionValue(CANCER_TRANS_FILE);
        SampleMutationsFile = cmd.getOptionValue(SAMPLE_MUT_FILE);

        Fusions = LoadTypes.contains(FUSION) ? new FusionCohortConfig(cmd) : null;

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
    }

    public static boolean isValid(final CommandLine cmd)
    {
        if(!cmd.hasOption(ROOT_DATA_DIRECTORY) || !cmd.hasOption(DATA_OUTPUT_DIR) || !cmd.hasOption(SAMPLE_DATA_FILE)
        || !cmd.hasOption(LOAD_TYPES))
        {
            return false;
        }

        return true;
    }

    public String formCohortFilename(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + "isofox_" + OutputIdentifier + "." + fileId;
        else
            return OutputDir + "isofox_" + fileId;
    }

    public static String formSampleFilename(final CohortConfig config, final String sampleId, final CohortAnalysisType dataType)
    {
        String filename = config.RootDataDir;

        if(config.UseSampleDirectories)
            filename += File.separator + sampleId + File.separator;

        filename += sampleId + ISOFOX_ID;
        filename += getFileId(dataType);
        return filename;
    }

    public static boolean formSampleFilenames(final CohortConfig config, final CohortAnalysisType dataType, final List<Path> filenames)
    {
        List<String> missingSampleIds = Lists.newArrayList();

        for(final String sampleId : config.SampleData.SampleIds)
        {
            String filename = formSampleFilename(config, sampleId, dataType);

            final Path path = Paths.get(filename);

            if (!Files.exists(path))
            {
                if(config.FailOnMissingSample)
                {
                    ISF_LOGGER.error("sampleId({}) file({}) not found", sampleId, filename);
                    filenames.clear();
                    return false;
                }
                else
                {
                    ISF_LOGGER.info("sampleId({}) file({}) not found, skipping", sampleId, filename);
                    missingSampleIds.add(sampleId);
                    continue;
                }
            }

            filenames.add(path);
        }

        config.SampleData.removeMissingSamples(missingSampleIds);

        return true;
    }

    public static Options createCmdLineOptions()
    {
        final Options options = new Options();
        options.addOption(ROOT_DATA_DIRECTORY, true, "Root data directory for input files or sample directories");
        options.addOption(SAMPLE_DATA_FILE, true, "File with list of samples and cancer types to load data for");
        options.addOption(DATA_OUTPUT_DIR, true, "Output directory");
        options.addOption(USE_SAMPLE_DIRS, false, "File with list of samples to load data for");
        options.addOption(FAIL_MISSING, false, "Exit if sample input file isn't found");
        options.addOption(ALL_AVAILABLE_FILES, false, "Load all files in root directory matching expeted Isofox file names");
        options.addOption(LOAD_TYPES, true, "List of data types to load & process");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(EXCLUDED_GENE_ID_FILE, true, "Optional CSV file of genes to ignore");
        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");


        options.addOption(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA, false, "Write per-sample gene distribution data file");
        options.addOption(CANCER_GENE_FILES, true, "Cancer gene distribution files, format: CancerType1-File1;CancerType2-File2");
        options.addOption(COHORT_TRANS_FILE, true, "Cohort transcript distribution file");
        options.addOption(CANCER_TRANS_FILE, true, "Cancer transcript distribution file");
        options.addOption(SAMPLE_MUT_FILE, true, "Sample mutations by gene and cancer type");
        options.addOption(TPM_ROUNDING, true, "TPM/FPM rounding factor, base-10 integer (default=2, ie 1%)");
        options.addOption(TPM_LOG_THRESHOLD, true, "Only write transcripts with TPM greater than this");

        AltSpliceJunctionCohort.addCmdLineOptions(options);
        FusionCohortConfig.addCmdLineOptions(options);

        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(THREADS, true, "Number of threads for task execution, default is 0 (off)");

        return options;
    }
}
