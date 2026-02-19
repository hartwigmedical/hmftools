package com.hartwig.hmftools.isofox.cohort;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.GENE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.EXPRESSION_DISTRIBUTION;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.EXTERNAL_EXPRESSION_COMPARE;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.FUSION;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.GENE_EXPRESSION_COMPARE;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.GENE_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.TRANSCRIPT_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.getIsofoxFileId;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig;
import com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortDistribution;
import com.hartwig.hmftools.isofox.expression.cohort.GeneratePanelNormalisation;
import com.hartwig.hmftools.isofox.fusion.cohort.FusionCohortConfig;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortAnalyser;
import com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortMatrix;
import com.hartwig.hmftools.isofox.novel.cohort.RecurrentVariantFinder;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache;
import com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher;
import com.hartwig.hmftools.isofox.unmapped.UmrCohortAnalyser;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class CohortConfig
{
    public static final String ROOT_DATA_DIRECTORY = "root_data_dir";
    public static final String SAMPLE_DATA_FILE = "sample_data_file";
    public static final String ALL_AVAILABLE_FILES = "all_available_files";
    public static final String ANALYSIS_TYPES = "analyses";
    public static final String FAIL_MISSING = "fail_on_missing_file";
    public static final String EXCLUDED_GENE_ID_FILE = "excluded_gene_id_file";

    public static final String SAMPLE_MUT_FILE = "sample_mut_file";

    public static final String COHORT_DELIM = CSV_DELIM;

    public final String RootDataDir;
    public final String OutputDir;
    public final String OutputIdentifier;
    public final SampleDataCache SampleData;
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

    public CohortConfig(final ConfigBuilder configBuilder)
    {
        RootDataDir = checkAddDirSeparator(configBuilder.getValue(ROOT_DATA_DIRECTORY, ""));;

        FailOnMissingSample = configBuilder.hasFlag(FAIL_MISSING);

        OutputDir = parseOutputDir(configBuilder);
        OutputIdentifier = configBuilder.getValue(OUTPUT_ID);

        final String sampleDataFile = configBuilder.getValue(SAMPLE_DATA_FILE);

        SampleData = new SampleDataCache(sampleDataFile);

        if(!SampleData.isValid())
        {
            ISF_LOGGER.warn("invalid sample data file({})", sampleDataFile);
        }

        AnalysisTypes = Arrays.stream(configBuilder.getValue(ANALYSIS_TYPES).split(ITEM_DELIM))
                .map(x -> AnalysisType.valueOf(x)).collect(Collectors.toList());

        RestrictedGeneIds = Lists.newArrayList();
        ExcludedGeneIds = Lists.newArrayList();

        if(configBuilder.hasValue(GENE_ID_FILE))
        {
            final String inputFile = configBuilder.getValue(GENE_ID_FILE);
            RestrictedGeneIds.addAll(loadGeneIdsFile(inputFile));
            ISF_LOGGER.info("file({}) loaded {} restricted genes", inputFile, RestrictedGeneIds.size());
        }

        if(configBuilder.hasValue(EXCLUDED_GENE_ID_FILE))
        {
            final String inputFile = configBuilder.getValue(EXCLUDED_GENE_ID_FILE);
            ExcludedGeneIds.addAll(loadGeneIdsFile(inputFile));
            ISF_LOGGER.info("file({}) loaded {} excluded genes", inputFile, ExcludedGeneIds.size());
        }

        EnsemblDataCache = configBuilder.getValue(ENSEMBL_DATA_DIR);

        RefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));

        ConvertUnmatchedCancerToOther = true;

        SampleMutationsFile = configBuilder.getValue(SAMPLE_MUT_FILE);

        Fusions = AnalysisTypes.contains(FUSION) ? new FusionCohortConfig(configBuilder) : null;

        Expression = requiresExpressionConfig(AnalysisTypes) ? new ExpressionCohortConfig(configBuilder) : null;

        DbAccess = createDatabaseAccess(configBuilder);

        Threads = parseThreads(configBuilder);
    }

    private static boolean requiresExpressionConfig(final List<AnalysisType> analysisTypes)
    {
        return analysisTypes.contains(GENE_EXPRESSION_COMPARE)
                || analysisTypes.contains(EXTERNAL_EXPRESSION_COMPARE) || analysisTypes.contains(GENE_EXPRESSION_MATRIX)
                || analysisTypes.contains(TRANSCRIPT_EXPRESSION_MATRIX) || analysisTypes.contains(EXPRESSION_DISTRIBUTION);
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
        String filename = convertWildcardSamplePath(config.RootDataDir, sampleId);
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

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(ROOT_DATA_DIRECTORY, true, "Root data directory for input files or sample directories");
        configBuilder.addConfigItem(SAMPLE_DATA_FILE, true, "File with list of samples and cancer types to load data for");
        configBuilder.addFlag(FAIL_MISSING, "Exit if sample input file isn't found");
        configBuilder.addConfigItem(ANALYSIS_TYPES, true, "List of data types to load & process");
        addEnsemblDir(configBuilder);
        configBuilder.addConfigItem(REF_GENOME, REF_GENOME_CFG_DESC);
        configBuilder.addPath(GENE_ID_FILE, false, GENE_ID_FILE_DESC);
        configBuilder.addPath(EXCLUDED_GENE_ID_FILE, false, "Optional CSV file of genes to ignore");
        configBuilder.addConfigItem(OUTPUT_ID, false, OUTPUT_ID_DESC);

        configBuilder.addConfigItem(SAMPLE_MUT_FILE, false, "Sample mutations by gene and cancer type");

        AltSjCohortAnalyser.registerConfig(configBuilder);
        SpliceVariantMatcher.registerConfig(configBuilder);
        FusionCohortConfig.registerConfig(configBuilder);
        ExpressionCohortConfig.registerConfig(configBuilder);
        ExpressionCohortDistribution.registerConfig(configBuilder);
        SpliceSiteCache.registerConfig(configBuilder);
        AltSjCohortMatrix.registerConfig(configBuilder);
        RecurrentVariantFinder.registerConfig(configBuilder);
        UmrCohortAnalyser.registerConfig(configBuilder);
        GeneratePanelNormalisation.registerConfig(configBuilder);

        addDatabaseCmdLineArgs(configBuilder, false);

        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}
