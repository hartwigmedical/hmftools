package com.hartwig.hmftools.lilac;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DATA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR_BAM_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FATAL_TOTAL_LOW_COVERAGE_POSITIONS;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FRAGS_PER_ALLELE;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FRAGS_REMOVE_SGL;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FREQUENCY_SCORE_PENALTY;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MAX_REF_FRAGMENTS;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_BASE_QUAL;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_DEPTH_FILTER;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE_SUPPORT;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_TOP_SCORE_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.LILAC_FILE_ID;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.ValidationStringency;

public class LilacConfig
{
    public final String Sample;
    public final String ReferenceBam;
    public final String TumorBam;
    public final String RnaBam;

    public final String ResourceDir;
    public final String RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final String SampleDataDir;
    public final String OutputDir;

    public final MhcClass ClassType;

    public final double HlaYPercentThreshold;

    public final int MinFragmentsPerAllele;
    public final int MinFragmentsToRemoveSingle;
    public final int MaxRefFragments;
    public final int Threads;
    public final double TopScoreThreshold;
    public final ValidationStringency BamStringency;

    public final String CopyNumberFile;
    public final String SomaticVariantsFile;

    public final boolean DebugPhasing;
    public final boolean RunValidation;
    public final int FatalTotalLowCoveragePositions;
    public final int MaxEliminationCandidates;
    public final boolean LogPerfCalcs;
    public final int StackSampleRate;

    // optional: pre-determine sample alleles, forced to be the final solution so coverage can be reported
    public final List<HlaAllele> ActualAlleles;

    // limit analysis from full data set to these + other configured alleles, for testing purposes
    public final List<HlaAllele> RestrictedAlleles;

    // config strings
    public static final String RESOURCE_DIR = "resource_dir";

    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String GENE_COPY_NUMBER = "gene_copy_number";

    public static final String MHC_CLASS = "mhc_class";

    // constant overrides
    private static final String MIN_BASE_QUAL = "min_base_qual";
    private static final String MAX_REF_FRAGMENTS = "max_ref_fragments";
    private static final String MIN_EVIDENCE_FACTOR = "min_evidence_factor";
    private static final String MIN_HIGH_QUAL_EVIDENCE_FACTOR = "min_high_qual_evidence_factor";
    private static final String MIN_EVIDENCE_SUPPORT = "min_evidence_support";
    private static final String MIN_DEPTH_FILTER = "min_depth_filter";
    private static final String MIN_FRAGMENTS_PER_ALLELE = "min_fragments_per_allele";
    private static final String MIN_FRAGMENTS_TO_REMOVE_SINGLE = "min_fragments_to_remove_single";
    private static final String TOP_SCORE_THRESHOLD = "top_score_threshold";
    private static final String HLA_Y_THRESHOLD = "hla_y_threshold";
    private static final String FREQUENCY_SCORE_PENALTY = "freq_score_penalty";

    // debug and technical config
    private static final String ACTUAL_ALLELES = "actual_alleles";
    private static final String RESTRICTED_ALLELES = "restricted_alleles";
    private static final String DEBUG_PHASING = "debug_phasing";
    private static final String RUN_VALIDATION = "run_validation";
    private static final String MAX_ELIM_CANDIDATES = "max_elim_candidates";
    private static final String FATAL_TOTAL_LOW_COVERAGE_POSITIONS = "fatal_total_low_coverage_positions";
    private static final String LOG_PERF_CALCS = "log_perf";
    private static final String STACK_SAMPLE_RATE = "stack_sample_rate";

    public static final Logger LL_LOGGER = LogManager.getLogger(LilacConfig.class);

    public LilacConfig(final ConfigBuilder configBuilder)
    {
        Sample = configBuilder.getValue(SAMPLE);

        if(configBuilder.hasValue(SAMPLE_DATA_DIR_CFG))
        {
            SampleDataDir = checkAddDirSeparator(configBuilder.getValue(SAMPLE_DATA_DIR_CFG));

            OutputDir = configBuilder.hasValue(OUTPUT_DIR) ? parseOutputDir(configBuilder) : SampleDataDir;

            String referenceId = Sample.substring(0, Sample.lastIndexOf('T')) + "R";
            ReferenceBam = SampleDataDir + referenceId + ".hla.bam";

            // these 3 are optional so check existence before initialising
            TumorBam = checkFileExists(SampleDataDir + Sample + ".hla.bam");
            RnaBam = checkFileExists(SampleDataDir + Sample + ".rna.hla.bam");
        }
        else
        {
            SampleDataDir = "";

            ReferenceBam = configBuilder.getValue(REFERENCE_BAM, "");
            TumorBam = configBuilder.getValue(TUMOR_BAM, "");
            RnaBam = configBuilder.getValue(RNA_BAM, "");

            if(configBuilder.hasValue(OUTPUT_DIR))
                OutputDir = parseOutputDir(configBuilder);
            else
                OutputDir = pathFromFile(ReferenceBam);
        }

        if(configBuilder.hasValue(SOMATIC_VCF) && configBuilder.hasValue(GENE_COPY_NUMBER))
        {
            SomaticVariantsFile = configBuilder.getValue(SOMATIC_VCF);
            CopyNumberFile = configBuilder.getValue(GENE_COPY_NUMBER);
        }
        else if(configBuilder.hasValue(PURPLE_DIR_CFG))
        {
            String purpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG));
            SomaticVariantsFile = PurpleCommon.purpleSomaticVcfFile(purpleDir, Sample);
            CopyNumberFile = GeneCopyNumberFile.generateFilename(purpleDir, Sample);
        }
        else if(!SampleDataDir.isEmpty())
        {
            String purpleGeneCopyNumberFile = GeneCopyNumberFile.generateFilename(SampleDataDir, Sample);
            CopyNumberFile = checkFileExists(purpleGeneCopyNumberFile);

            String somaticVariantsFile = PurpleCommon.purpleSomaticVcfFile(SampleDataDir, Sample);
            SomaticVariantsFile = checkFileExists(somaticVariantsFile);
        }
        else
        {
            SomaticVariantsFile = "";
            CopyNumberFile = "";
        }

        ResourceDir = checkAddDirSeparator(configBuilder.getValue(RESOURCE_DIR));
        RefGenome = configBuilder.getValue(REF_GENOME, "");

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        ClassType = MhcClass.valueOf(configBuilder.getValue(MHC_CLASS));

        LilacConstants.LOW_BASE_QUAL_THRESHOLD = (byte)configBuilder.getInteger(MIN_BASE_QUAL);
        LilacConstants.MIN_EVIDENCE_FACTOR = configBuilder.getDecimal(MIN_EVIDENCE_FACTOR);
        LilacConstants.MIN_HIGH_QUAL_EVIDENCE_FACTOR = configBuilder.getDecimal(MIN_HIGH_QUAL_EVIDENCE_FACTOR);
        LilacConstants.MIN_EVIDENCE_SUPPORT = configBuilder.getInteger(MIN_EVIDENCE_SUPPORT);
        LilacConstants.MIN_DEPTH_FILTER = configBuilder.getInteger(MIN_DEPTH_FILTER);

        LilacConstants.FREQUENCY_SCORE_PENALTY = configBuilder.getDecimal(FREQUENCY_SCORE_PENALTY);

        MaxRefFragments = configBuilder.getInteger(MAX_REF_FRAGMENTS);
        HlaYPercentThreshold = configBuilder.getDecimal(HLA_Y_THRESHOLD);

        MinFragmentsPerAllele = configBuilder.getInteger(MIN_FRAGMENTS_PER_ALLELE);
        MinFragmentsToRemoveSingle = configBuilder.getInteger(MIN_FRAGMENTS_TO_REMOVE_SINGLE);
        FatalTotalLowCoveragePositions = configBuilder.getInteger(FATAL_TOTAL_LOW_COVERAGE_POSITIONS);

        TopScoreThreshold = min(configBuilder.getDecimal(TOP_SCORE_THRESHOLD), 0.5);

        ActualAlleles = parseAlleleList(configBuilder.getValue(ACTUAL_ALLELES));
        RestrictedAlleles = parseAlleleList(configBuilder.getValue(RESTRICTED_ALLELES));
        MaxEliminationCandidates = configBuilder.getInteger(MAX_ELIM_CANDIDATES);

        Threads = parseThreads(configBuilder);
        BamStringency = BamUtils.validationStringency(configBuilder);

        DebugPhasing = configBuilder.hasFlag(DEBUG_PHASING);
        RunValidation = configBuilder.hasFlag(RUN_VALIDATION);
        LogPerfCalcs = configBuilder.hasFlag(LOG_PERF_CALCS);
        StackSampleRate = configBuilder.getInteger(STACK_SAMPLE_RATE);

        if(!checkCreateOutputDir(OutputDir))
        {
            LL_LOGGER.error("failed to create output directory: {}", OutputDir);
            System.exit(1);
        }
    }

    public boolean tumorOnly() { return ReferenceBam.isEmpty() && !TumorBam.isEmpty(); }
    public String alleleCallingBam() { return !ReferenceBam.isEmpty() ? ReferenceBam : TumorBam; }

    private static String checkFileExists(final String filename)
    {
        return Files.exists(Paths.get(filename)) ? filename : "";
    }

    public String formFileId(final String fileId) { return OutputDir + Sample + LILAC_FILE_ID + fileId; }

    public void logParams()
    {
        LL_LOGGER.info("sample({}) inputs: referenceBam({}) tumorBam({}) somaticVCF({}) geneCopyNumber({}) rnaBam({})",
                Sample, !ReferenceBam.isEmpty(), !TumorBam.isEmpty(), !SomaticVariantsFile.isEmpty(), !CopyNumberFile.isEmpty(), !RnaBam.isEmpty());

        if(RunValidation)
        {
            LL_LOGGER.info("running validation routines");
        }

        if(!ActualAlleles.isEmpty())
        {
            LL_LOGGER.info("expected alleles: {}", HlaAllele.toString(ActualAlleles));
        }

        if(!RestrictedAlleles.isEmpty())
        {
            LL_LOGGER.info("restricted alleles: {}", HlaAllele.toString(RestrictedAlleles));
        }
    }

    public LilacConfig(final String sampleId)
    {
        OutputDir = "";
        Sample = sampleId;
        ReferenceBam = "";
        TumorBam = "";
        RnaBam = "";
        ResourceDir = "";
        SampleDataDir = "";
        RefGenome = "";
        RefGenVersion = V37;

        ClassType = MhcClass.CLASS_1;

        MaxRefFragments = DEFAULT_MAX_REF_FRAGMENTS;

        MinFragmentsPerAllele = DEFAULT_FRAGS_PER_ALLELE;
        MinFragmentsToRemoveSingle = DEFAULT_FRAGS_REMOVE_SGL;
        TopScoreThreshold = DEFAULT_TOP_SCORE_THRESHOLD;
        FatalTotalLowCoveragePositions = DEFAULT_FATAL_TOTAL_LOW_COVERAGE_POSITIONS * 2;
        HlaYPercentThreshold = DEFAULT_HLA_Y_FRAGMENT_THRESHOLD;

        CopyNumberFile = "";
        SomaticVariantsFile = "";

        ActualAlleles = Lists.newArrayList();
        RestrictedAlleles = Lists.newArrayList();
        MaxEliminationCandidates = 0;

        BamStringency = ValidationStringency.STRICT;
        Threads = 0;
        DebugPhasing = false;
        RunValidation = true;
        LogPerfCalcs = false;
        StackSampleRate = 0;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addRequiredConfigItem(SAMPLE, "Name of sample");
        configBuilder.addPath(SAMPLE_DATA_DIR_CFG, false, SAMPLE_DATA_DIR_DESC);

        configBuilder.addPath(REFERENCE_BAM, false, REFERENCE_BAM_DESC);
        configBuilder.addPath(TUMOR_BAM, false, TUMOR_BAM_DESC);
        configBuilder.addPath(RNA_BAM, false, RNA_BAM_DESC);

        registerCommonConfig(configBuilder);

        configBuilder.addInteger(MIN_BASE_QUAL, "Min base quality threshold", DEFAULT_MIN_BASE_QUAL);
        configBuilder.addDecimal(MIN_EVIDENCE_FACTOR, "Min fragment evidence required", DEFAULT_MIN_EVIDENCE_FACTOR);
        configBuilder.addInteger(MAX_REF_FRAGMENTS, "Cap ref fragments in solution search, 0 uses all", DEFAULT_MAX_REF_FRAGMENTS);
        configBuilder.addDecimal(MIN_HIGH_QUAL_EVIDENCE_FACTOR, "Min high-qual fragment evidence factor", DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR);
        configBuilder.addInteger(MIN_EVIDENCE_SUPPORT, "Min fragment evidence depth", DEFAULT_MIN_EVIDENCE_SUPPORT);
        configBuilder.addInteger(MIN_DEPTH_FILTER, "Raw depth below which we do not filter", DEFAULT_MIN_DEPTH_FILTER);
        configBuilder.addInteger(MIN_FRAGMENTS_PER_ALLELE, "Min fragments per allele", DEFAULT_FRAGS_PER_ALLELE);
        configBuilder.addInteger(MIN_FRAGMENTS_TO_REMOVE_SINGLE, "Min fragments to remote single", DEFAULT_FRAGS_REMOVE_SGL);

        configBuilder.addInteger(FATAL_TOTAL_LOW_COVERAGE_POSITIONS, "Fatal total low coverage positions across all HLA loci", DEFAULT_FATAL_TOTAL_LOW_COVERAGE_POSITIONS);
        configBuilder.addDecimal(HLA_Y_THRESHOLD, "HLA-Y percent threshold", DEFAULT_HLA_Y_FRAGMENT_THRESHOLD);
        configBuilder.addDecimal(FREQUENCY_SCORE_PENALTY, "Population frequency penalty", DEFAULT_FREQUENCY_SCORE_PENALTY);

        configBuilder.addDecimal(TOP_SCORE_THRESHOLD, "Max distance from top score", DEFAULT_TOP_SCORE_THRESHOLD);
        configBuilder.addConfigItem(ACTUAL_ALLELES, "Comma separated known actual alleles for the sample");
        configBuilder.addConfigItem(RESTRICTED_ALLELES, "Comma separated restricted analysis allele list");

        configBuilder.addInteger(
                MAX_ELIM_CANDIDATES,
                "Revert to only common alleles if candidate allele count exceeds this after elimination", 0);

        configBuilder.addPath(GENE_COPY_NUMBER, false, "Path to gene copy number file");
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(SOMATIC_VCF, false, "Path to sample Purple somatic VCF");
        configBuilder.addFlag(DEBUG_PHASING, "More detailed logging of phasing");
        configBuilder.addFlag(RUN_VALIDATION, "Run validation checks");
        configBuilder.addFlag(LOG_PERF_CALCS, "Log performance metrics");
        configBuilder.addInteger(STACK_SAMPLE_RATE, "Stack sampling rate", 0);
        ResultsWriter.registerConfig(configBuilder);

        BamUtils.addValidationStringencyOption(configBuilder);

        addRefGenomeConfig(configBuilder, true);
        addOutputDir(configBuilder);
        addThreadOptions(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public static void registerCommonConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(RESOURCE_DIR, true, "Path to resource files");

        configBuilder.addConfigItem(MHC_CLASS, false, "MHC Class Type", MhcClass.CLASS_1.toString());
    }

    private static List<HlaAllele> parseAlleleList(final String allelesStr)
    {
        if(allelesStr == null || allelesStr.isEmpty())
            return Lists.newArrayList();

        String[] alleles = allelesStr.split(ITEM_DELIM, -1);
        return Arrays.stream(alleles).map(HlaAllele::fromString).collect(Collectors.toList());
    }
}
