package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.isofox.IsofoxConfig.DATA_OUTPUT_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.EXCLUDED_GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_ID_FILE;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.OUTPUT_ID;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.EXPRESSION_DISTRIBUTION;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.EXTERNAL_EXPRESSION_COMPARE;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.GENE_EXPRESSION_COMPARE;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.GENE_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.TRANSCRIPT_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.getIsofoxFileId;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortDistribution;
import com.hartwig.hmftools.isofox.fusion.cohort.FusionCohortConfig;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortAnalyser;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortMatrix;
import com.hartwig.hmftools.isofox.novel.cohort.RecurrentVariantFinder;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class CohortConfig
{
    public static final String ROOT_DATA_DIRECTORY = "root_data_dir";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    public static final String USE_SAMPLE_DIRS = "use_sample_dir";
    public static final String ALL_AVAILABLE_FILES = "all_available_files";
    public static final String ANALYSIS_TYPES = "analyses";
    public static final String FAIL_MISSING = "fail_on_missing_file";

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

    public final List<AnalysisType> AnalysisTypes;

    public final DatabaseAccess DbAccess;

    public final String EnsemblDataCache;
    public final RefGenomeInterface RefGenome;

    public final String SampleMutationsFile;

    public final FusionCohortConfig Fusions;

    public final ExpressionCohortConfig Expression;

    public final int Threads;

    public CohortConfig(final CommandLine cmd)
    {
        RootDataDir = cmd.hasOption(ROOT_DATA_DIRECTORY) ? checkAddDirSeparator(cmd.getOptionValue(ROOT_DATA_DIRECTORY)) : "";;

        UseSampleDirectories = cmd.hasOption(USE_SAMPLE_DIRS);
        AllAvailableFiles = !UseSampleDirectories && cmd.hasOption(ALL_AVAILABLE_FILES);
        FailOnMissingSample = cmd.hasOption(FAIL_MISSING);

        OutputDir = parseOutputDir(cmd);
        OutputIdentifier = cmd.getOptionValue(OUTPUT_ID);

        final String sampleDataFile = cmd.getOptionValue(SAMPLE_DATA_FILE);

        SampleData = new SampleDataCache(sampleDataFile);

        if(!SampleData.isValid())
        {
            ISF_LOGGER.warn("invalid sample data file({})", sampleDataFile);
        }

        AnalysisTypes = Arrays.stream(cmd.getOptionValue(ANALYSIS_TYPES).split(ITEM_DELIM))
                .map(x -> AnalysisType.valueOf(x)).collect(Collectors.toList());

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));
            ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }

        if(cmd.hasOption(EXCLUDED_GENE_ID_FILE))
        {
            final String inputFile = cmd.getOptionValue(EXCLUDED_GENE_ID_FILE);
            ExcludedGeneIds.addAll(loadGeneIdsFile(inputFile));
            ISF_LOGGER.info("file({}) loaded {} excluded genes", inputFile, ExcludedGeneIds.size());
        }

        EnsemblDataCache = cmd.getOptionValue(GENE_TRANSCRIPTS_DIR);

        RefGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        ConvertUnmatchedCancerToOther = true;

        SampleMutationsFile = cmd.getOptionValue(SAMPLE_MUT_FILE);

        Fusions = AnalysisTypes.contains(FUSION) ? new FusionCohortConfig(cmd) : null;

        Expression = requiresExpressionConfig(AnalysisTypes) ? new ExpressionCohortConfig(cmd) : null;

        DbAccess = createDatabaseAccess(cmd);

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
    }

    private static boolean requiresExpressionConfig(final List<AnalysisType> analysisTypes)
    {
        return analysisTypes.contains(GENE_EXPRESSION_COMPARE)
                || analysisTypes.contains(EXTERNAL_EXPRESSION_COMPARE) || analysisTypes.contains(GENE_EXPRESSION_MATRIX)
                || analysisTypes.contains(TRANSCRIPT_EXPRESSION_MATRIX) || analysisTypes.contains(EXPRESSION_DISTRIBUTION);
    }

    public static boolean isValid(final CommandLine cmd)
    {
        if(!cmd.hasOption(DATA_OUTPUT_DIR) || !cmd.hasOption(SAMPLE_DATA_FILE) || !cmd.hasOption(ANALYSIS_TYPES))
        {
            return false;
        }

        return true;
    }

    public String formCohortFilename(final String fileId)
    {
        if(OutputIdentifier != null)
            return OutputDir + "isofox." + OutputIdentifier + "." + fileId;
        else
            return OutputDir + "isofox." + fileId;
    }

    public static String formSampleFilename(final CohortConfig config, final String sampleId, final AnalysisType dataType)
    {
        String filename = config.RootDataDir;

        if(config.UseSampleDirectories)
            filename += File.separator + sampleId + File.separator;

        filename += sampleId + ISF_FILE_ID;
        filename += getIsofoxFileId(dataType);
        return filename;
    }

    public static boolean formSampleFilenames(final CohortConfig config, final AnalysisType dataType, final List<Path> filenames)
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
        options.addOption(ANALYSIS_TYPES, true, "List of data types to load & process");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Path to Ensembl data cache");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(GENE_ID_FILE, true, "Optional CSV file of genes to analyse");
        options.addOption(EXCLUDED_GENE_ID_FILE, true, "Optional CSV file of genes to ignore");
        options.addOption(OUTPUT_ID, true, "Optionally add identifier to output files");

        options.addOption(SAMPLE_MUT_FILE, true, "Sample mutations by gene and cancer type");

        AltSjCohortAnalyser.addCmdLineOptions(options);
        SpliceVariantMatcher.addCmdLineOptions(options);
        FusionCohortConfig.addCmdLineOptions(options);
        ExpressionCohortConfig.addCmdLineOptions(options);
        ExpressionCohortDistribution.addCmdLineOptions(options);
        SpliceSiteCache.addCmdLineOptions(options);
        AltSjCohortMatrix.addCmdLineOptions(options);
        RecurrentVariantFinder.addCmdLineOptions(options);

        addDatabaseCmdLineArgs(options);

        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(THREADS, true, "Number of threads for task execution, default is 0 (off)");

        return options;
    }
}
