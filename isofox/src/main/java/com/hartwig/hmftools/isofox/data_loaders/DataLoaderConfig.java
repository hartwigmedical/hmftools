package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.isofox.IsofoxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.EXCLUDED_GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.LOG_DEBUG;
import static com.hartwig.hmftools.isofox.IsofoxConfig.OUTPUT_ID;
import static com.hartwig.hmftools.isofox.IsofoxConfig.loadGeneIdsFile;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoadType.getFileId;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import jdk.nashorn.internal.ir.EmptyNode;

public class DataLoaderConfig
{
    public static final String ROOT_DATA_DIRECTORY = "root_data_dir";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    public static final String USE_SAMPLE_DIRS = "use_sample_dir";
    public static final String ALL_AVAILABLE_FILES = "all_available_files";
    public static final String LOAD_TYPES = "load_types";
    public static final String ALT_SJ_MIN_SAMPLES = "alt_sj_min_samples";
    public static final String ALT_SJ_PROB_THRESHOLD = "alt_sj_prob_threshold";
    public static final String ALT_SJ_MIN_FRAGS_REQ_GENES = "alt_sj_min_frags_req_genes";
    public static final String SPLICE_VARIANT_FILE = "splice_variant_file";

    public static final String TPM_LOG_THRESHOLD = "tpm_log_threshold";
    public static final String TPM_ROUNDING = "tpm_rounding";

    public static final String WRITE_SAMPLE_GENE_DISTRIBUTION_DATA = "write_sample_gene_dist";

    public static final String CANCER_GENE_FILES = "cancer_gene_files";

    public static final String COHORT_TRANS_FILE = "cohort_trans_file";
    public static final String CANCER_TRANS_FILE = "cancer_trans_file";

    public final String RootDataDir;
    public final String OutputDir;
    public final String OutputIdentifier;
    public final SampleDataCache SampleData;
    public final boolean UseSampleDirectories;
    public final boolean AllAvailableFiles;
    public final List<String> RestrictedGeneIds;
    public final List<String> ExcludedGeneIds;

    public final List<DataLoadType> LoadTypes;

    public final int AltSJMinSampleThreshold;
    public final int AltSJMinFragsUnrequiredGenes;
    public final double AltSJProbabilityThreshold;
    public final String SpliceVariantFile;
    public final String EnsemblDataCache;

    public final boolean WriteSampleGeneDistributionData;

    public final String CohortTransFile;
    public final String CancerTransFile;

    public final String CancerGeneFiles;
    public final double TpmLogThreshold;
    public final double TpmRounding;

    public DataLoaderConfig(final CommandLine cmd)
    {
        RootDataDir = cmd.getOptionValue(ROOT_DATA_DIRECTORY);
        UseSampleDirectories = cmd.hasOption(USE_SAMPLE_DIRS);
        AllAvailableFiles = !UseSampleDirectories && cmd.hasOption(ALL_AVAILABLE_FILES);

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

        AltSJMinSampleThreshold = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_SAMPLES, "0"));
        AltSJMinFragsUnrequiredGenes = Integer.parseInt(cmd.getOptionValue(ALT_SJ_MIN_FRAGS_REQ_GENES, "0"));
        AltSJProbabilityThreshold = Double.parseDouble(cmd.getOptionValue(ALT_SJ_PROB_THRESHOLD, "1.0"));
        SpliceVariantFile = cmd.getOptionValue(SPLICE_VARIANT_FILE);
        EnsemblDataCache = cmd.getOptionValue(GENE_TRANSCRIPTS_DIR);

        TpmLogThreshold = Double.parseDouble(cmd.getOptionValue(TPM_LOG_THRESHOLD, "0"));
        TpmRounding = Double.parseDouble(cmd.getOptionValue(TPM_ROUNDING, "2"));

        WriteSampleGeneDistributionData = cmd.hasOption(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA);
        CohortTransFile = cmd.getOptionValue(COHORT_TRANS_FILE);
        CancerTransFile = cmd.getOptionValue(CANCER_TRANS_FILE);

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
    }

    public static boolean isValid(final CommandLine cmd)
    {
        return cmd.hasOption(ROOT_DATA_DIRECTORY) && cmd.hasOption(DATA_OUTPUT_DIR) && cmd.hasOption(SAMPLE_DATA_FILE)
                && cmd.hasOption(LOAD_TYPES);
    }

    public String formCohortFilename(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + "isofox_" + OutputIdentifier + "." + fileId;
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
        options.addOption(ALL_AVAILABLE_FILES, false, "Load all files in root directory matching expeted Isofox file names");
        options.addOption(LOAD_TYPES, true, "List of data types to load & process");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache");
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(EXCLUDED_GENE_ID_FILE, true, "Optional CSV file of genes to ignore");
        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");

        options.addOption(ALT_SJ_MIN_SAMPLES, true, "Min number of samples to report an alt SJ");
        options.addOption(ALT_SJ_MIN_FRAGS_REQ_GENES, true, "Min frag count supporting alt-SJs outside gene panel");
        options.addOption(ALT_SJ_PROB_THRESHOLD, true, "Only write alt SJs for fisher probability less than this");
        options.addOption(SPLICE_VARIANT_FILE, true, "File with somatic variants potentially affecting splicing");

        options.addOption(WRITE_SAMPLE_GENE_DISTRIBUTION_DATA, false, "Write per-sample gene distribution data file");
        options.addOption(CANCER_GENE_FILES, true, "Cancer gene distribution files, format: CancerType1-File1;CancerType2-File2");
        options.addOption(COHORT_TRANS_FILE, true, "Cohort transcript distribution file");
        options.addOption(CANCER_TRANS_FILE, true, "Cancer transcript distribution file");
        options.addOption(TPM_ROUNDING, true, "TPM/FPM rounding factor, base-10 integer (default=2, ie 1%)");
        options.addOption(TPM_LOG_THRESHOLD, true, "Only write transcripts with TPM greater than this");

        options.addOption(LOG_DEBUG, false, "Log verbose");

        return options;
    }
}
