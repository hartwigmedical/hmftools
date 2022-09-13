package com.hartwig.hmftools.lilac;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeConfig;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FRAGS_PER_ALLELE;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FRAGS_REMOVE_SGL;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_BASE_QUAL;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_TOP_SCORE_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FATAL_LOW_COVERAGE_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_HLA_Y_FRAGMENT_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.ITEM_DELIM;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LilacConfig
{
    private final CommandLine mCmdLineArgs;

    public final String Sample;
    public final String ReferenceBam;
    public final String TumorBam;
    public final String RnaBam;
    public final String RunId;

    public final String ResourceDir;
    public final String RefGenome;
    public final RefGenomeVersion RefGenVersion;
    public final String SampleDataDir;
    public final String OutputDir;
    public final boolean WriteAllFiles;

    public int MinBaseQual;
    private final int MinEvidence;
    private final double MinEvidenceFactor;
    private final double MinHighQualEvidenceFactor;
    public final double HlaYPercentThreshold;

    public final int MinFragmentsPerAllele;
    public final int MinFragmentsToRemoveSingle;
    public final int Threads;
    public final double TopScoreThreshold;

    public final String CopyNumberFile;
    public final String SomaticVariantsFile;
    public final String CohortSomaticVariantsFile;

    public final boolean DebugPhasing;
    public final boolean RunValidation;
    public final int FatalLowCoverage;
    public final int MaxEliminationCandidates;
    public final boolean LogPerfCalcs;

    // optional: pre-determine sample alleles, forced to be the final solution so coverage can be reported
    public final List<HlaAllele> ActualAlleles;

    // limit analysis from full data set to these + other configured alleles, for testing purposes
    public final List<HlaAllele> RestrictedAlleles;

    // config strings
    public static final String RESOURCE_DIR = "resource_dir";

    private static final String SAMPLE = "sample";
    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String RNA_BAM = "rna_bam";
    private static final String RUN_ID = "run_id";

    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String COHORT_SOMATIC_VCF = "cohort_somatic_file";
    private static final String GENE_COPY_NUMBER = "gene_copy_number";

    // constant overrides
    private static final String MIN_BASE_QUAL = "min_base_qual";
    private static final String MIN_EVIDENCE = "min_evidence";
    private static final String MIN_EVIDENCE_FACTOR = "min_evidence_factor";
    private static final String MIN_HIGH_QUAL_EVIDENCE_FACTOR = "min_high_qual_evidence_factor";
    private static final String MIN_FRAGMENTS_PER_ALLELE = "min_fragments_per_allele";
    private static final String MIN_FRAGMENTS_TO_REMOVE_SINGLE = "min_fragments_to_remove_single";
    private static final String TOP_SCORE_THRESHOLD = "top_score_threshold";
    private static final String HLA_Y_THRESHOLD = "hla_y_threshold";

    // debug and technical config
    private static final String ACTUAL_ALLELES = "actual_alleles";
    private static final String RESTRICTED_ALLELES = "restricted_alleles";
    private static final String DEBUG_PHASING = "debug_phasing";
    public static final String LOG_DEBUG = "log_debug";
    public static final String RUN_VALIDATION = "run_validation";
    public static final String MAX_ELIM_CANDIDATES = "max_elim_candidates";
    public static final String WRITE_ALL_FILES = "write_all_files";
    public static final String FATAL_LOW_COVERAGE = "fatal_low_coverage";
    public static final String LOG_PERF_CALCS = "log_perf";

    public static final Logger LL_LOGGER = LogManager.getLogger(LilacConfig.class);;

    public LilacConfig(final CommandLine cmd)
    {
        mCmdLineArgs = cmd;

        Sample = cmd.getOptionValue(SAMPLE);

        if(cmd.hasOption(SAMPLE_DATA_DIR))
        {
            SampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR));

            OutputDir = SampleDataDir;

            String referenceId = Sample.substring(0, Sample.lastIndexOf('T')) + "R";
            ReferenceBam = SampleDataDir + referenceId + ".hla.bam";

            // these 3 are optional so check existence before initialising
            TumorBam = checkFileExists(SampleDataDir + Sample + ".hla.bam");
            RnaBam = checkFileExists(SampleDataDir + Sample + ".rna.hla.bam");

            String purpleGeneCopyNumberFile = GeneCopyNumberFile.generateFilenameForReading(SampleDataDir, Sample);
            CopyNumberFile = checkFileExists(purpleGeneCopyNumberFile);

            String somaticVariantsFile = PurpleCommon.purpleSomaticVcfFile(SampleDataDir, Sample);
            SomaticVariantsFile = checkFileExists(somaticVariantsFile);
            CohortSomaticVariantsFile = "";
        }
        else
        {
            SampleDataDir = "";

            if(!cmd.hasOption(REFERENCE_BAM) && cmd.hasOption(TUMOR_BAM))
            {
                // interpret this as tumor-only mode
                ReferenceBam = cmd.getOptionValue(TUMOR_BAM, "");
                TumorBam = "";
            }
            else
            {
                ReferenceBam = cmd.getOptionValue(REFERENCE_BAM, "");
                TumorBam = cmd.getOptionValue(TUMOR_BAM, "");
            }

            RnaBam = cmd.getOptionValue(RNA_BAM, "");

            SomaticVariantsFile = cmd.getOptionValue(SOMATIC_VCF, "");
            CohortSomaticVariantsFile = cmd.getOptionValue(COHORT_SOMATIC_VCF, "");
            CopyNumberFile = cmd.getOptionValue(GENE_COPY_NUMBER, "");

            OutputDir = parseOutputDir(cmd);
        }

        RunId = cmd.getOptionValue(RUN_ID, "");

        ResourceDir = checkAddDirSeparator(cmd.getOptionValue(RESOURCE_DIR));
        RefGenome = cmd.getOptionValue(REF_GENOME, "");

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;

        MinBaseQual = getConfigValue(cmd, MIN_BASE_QUAL, DEFAULT_MIN_BASE_QUAL);
        MinEvidence = getConfigValue(cmd, MIN_EVIDENCE, DEFAULT_MIN_EVIDENCE);
        MinEvidenceFactor = getConfigValue(cmd, MIN_EVIDENCE_FACTOR, DEFAULT_MIN_EVIDENCE_FACTOR);
        MinHighQualEvidenceFactor = getConfigValue(cmd, MIN_HIGH_QUAL_EVIDENCE_FACTOR, DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR);
        HlaYPercentThreshold = getConfigValue(cmd, HLA_Y_THRESHOLD, DEFAULT_HLA_Y_FRAGMENT_THRESHOLD);

        MinFragmentsPerAllele = getConfigValue(cmd, MIN_FRAGMENTS_PER_ALLELE, DEFAULT_FRAGS_PER_ALLELE);
        MinFragmentsToRemoveSingle = getConfigValue(cmd, MIN_FRAGMENTS_TO_REMOVE_SINGLE, DEFAULT_FRAGS_REMOVE_SGL);
        FatalLowCoverage = getConfigValue(cmd, FATAL_LOW_COVERAGE, DEFAULT_FATAL_LOW_COVERAGE_THRESHOLD);

        TopScoreThreshold = min(getConfigValue(cmd, TOP_SCORE_THRESHOLD, DEFAULT_TOP_SCORE_THRESHOLD), 0.5);

        ActualAlleles = parseAlleleList(cmd.getOptionValue(ACTUAL_ALLELES));
        RestrictedAlleles = parseAlleleList(cmd.getOptionValue(RESTRICTED_ALLELES));
        MaxEliminationCandidates = Integer.parseInt(cmd.getOptionValue(MAX_ELIM_CANDIDATES, "0"));

        Threads = parseThreads(cmd);

        DebugPhasing = cmd.hasOption(DEBUG_PHASING);
        RunValidation = cmd.hasOption(RUN_VALIDATION);
        WriteAllFiles = cmd.hasOption(WRITE_ALL_FILES);
        LogPerfCalcs = cmd.hasOption(LOG_PERF_CALCS);
    }

    private String checkFileExists(final String filename)
    {
        return Files.exists(Paths.get(filename)) ? filename : "";
    }

    public String formFileId(final String fileId) { return OutputDir + Sample + "." + fileId; }

    public double calcMinEvidence(int totalFragments) { return max(MinEvidence, totalFragments * MinEvidenceFactor); }
    public double calcMinHighQualEvidence(int totalFragments) { return totalFragments * MinHighQualEvidenceFactor; }

    public boolean isValid()
    {
        if(mCmdLineArgs == null)
            return true;

        if(ReferenceBam.isEmpty() || !Files.exists(Paths.get(ReferenceBam)))
        {
            LL_LOGGER.error("missing or invalid reference BAM");
            return false;
        }

        if(RefGenome.isEmpty() || !Files.exists(Paths.get(ResourceDir)))
        {
            LL_LOGGER.error("missing or invalid reference genome");
            return false;
        }

        if(ResourceDir.isEmpty() || !Files.exists(Paths.get(ResourceDir)))
        {
            LL_LOGGER.error("missing resource file directory");
            return false;
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            LL_LOGGER.error("failed to create output directory: {}", OutputDir);
            return false;
        }

        return true;
    }

    public void logParams()
    {
        LL_LOGGER.info("sample({}) inputs: tumorBam({}) somaticVCF({}) geneCopyNumber({}) rnaBam({})",
                Sample, !TumorBam.isEmpty(), !SomaticVariantsFile.isEmpty(), !CopyNumberFile.isEmpty(), !RnaBam.isEmpty());

        LL_LOGGER.info("minBaseQual({}), minEvidence({}) minFragmentsPerAllele({}) "
                + "minFragmentsToRemoveSingle({}) maxDistanceFromTopScore({})",
                MinBaseQual, MinEvidence, MinFragmentsPerAllele, MinFragmentsToRemoveSingle, TopScoreThreshold);

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
        mCmdLineArgs = null;

        OutputDir = "";
        Sample = sampleId;
        ReferenceBam = "";
        TumorBam = "";
        RnaBam = "";
        ResourceDir = "";
        SampleDataDir = "";
        RefGenome = "";
        RefGenVersion = V37;
        RunId = "";

        MinBaseQual = DEFAULT_MIN_BASE_QUAL;
        MinEvidence = DEFAULT_MIN_EVIDENCE;
        MinEvidenceFactor = DEFAULT_MIN_EVIDENCE_FACTOR;
        MinHighQualEvidenceFactor = DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR;

        MinFragmentsPerAllele = DEFAULT_FRAGS_PER_ALLELE;
        MinFragmentsToRemoveSingle = DEFAULT_FRAGS_REMOVE_SGL;
        TopScoreThreshold = DEFAULT_TOP_SCORE_THRESHOLD;
        FatalLowCoverage = DEFAULT_FATAL_LOW_COVERAGE_THRESHOLD;
        HlaYPercentThreshold = DEFAULT_HLA_Y_FRAGMENT_THRESHOLD;

        CopyNumberFile = "";
        SomaticVariantsFile = "";
        CohortSomaticVariantsFile = "";

        ActualAlleles = Lists.newArrayList();
        RestrictedAlleles = Lists.newArrayList();
        MaxEliminationCandidates = 0;

        Threads = 0;
        DebugPhasing = false;
        RunValidation = true;
        WriteAllFiles = false;
        LogPerfCalcs = false;
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(SAMPLE_DATA_DIR, true,"Path to all sample files");
        options.addOption(REFERENCE_BAM, true,"Path to reference/normal BAM");
        options.addOption(TUMOR_BAM, true,"Path to tumor BAM");
        options.addOption(RNA_BAM, true,"Analyse tumor BAM only");
        options.addOption(RUN_ID, true,"Only search for HLA-Y fragments");
        options.addOption(RESOURCE_DIR, true,"Path to resource files");

        options.addOption(MIN_BASE_QUAL, true,"Min base quality threshold, default: " + DEFAULT_MIN_BASE_QUAL);

        options.addOption(MIN_EVIDENCE, true,"Min fragment evidence required, default: " + DEFAULT_MIN_EVIDENCE);

        options.addOption(
                MIN_HIGH_QUAL_EVIDENCE_FACTOR, true,
                "Min high-qual fragment evidence factor, default: " + DEFAULT_MIN_HIGH_QUAL_EVIDENCE_FACTOR);

        options.addOption(
                MIN_EVIDENCE_FACTOR, true,"Min fragment evidence factor, default: " + DEFAULT_MIN_EVIDENCE_FACTOR);

        options.addOption(
                MIN_FRAGMENTS_PER_ALLELE, true,"Min fragments per allele, default: " + DEFAULT_FRAGS_PER_ALLELE);

        options.addOption(
                MIN_FRAGMENTS_TO_REMOVE_SINGLE, true,"Min fragments to remote single, default: " + DEFAULT_FRAGS_REMOVE_SGL);

        options.addOption(
                FATAL_LOW_COVERAGE, true,"Fatal low coverage, default: " + DEFAULT_FATAL_LOW_COVERAGE_THRESHOLD);

        options.addOption(
                HLA_Y_THRESHOLD, true,"HLA-Y percent threshold, default: " + DEFAULT_HLA_Y_FRAGMENT_THRESHOLD);

        options.addOption(TOP_SCORE_THRESHOLD, true,"Max distance from top score");
        options.addOption(ACTUAL_ALLELES, true,"Comma separated known actual alleles for the sample");
        options.addOption(RESTRICTED_ALLELES, true,"Comma separated restricted analysis allele list");
        options.addOption(MAX_ELIM_CANDIDATES, true, "Revert to only common alleles if candidate allele count exceeds this after elimination");
        options.addOption(GENE_COPY_NUMBER, true,"Deprecated, use 'gene_copy_number instead': Path to gene copy number file ");
        options.addOption(SOMATIC_VCF, true,"Path to somatic VCF");
        options.addOption(COHORT_SOMATIC_VCF, true,"Cohort somatic variants file");
        options.addOption(DEBUG_PHASING, false, "More detailed logging of phasing");
        options.addOption(RUN_VALIDATION, false, "Run validation checks");
        options.addOption(WRITE_ALL_FILES, false,"Write more detailed output files");
        options.addOption(LOG_PERF_CALCS, false,"Log performance metrics");

        addRefGenomeConfig(options);
        addOutputDir(options);
        addLoggingOptions(options);
        addThreadOptions(options);
        return options;
    }

    private final List<HlaAllele> parseAlleleList(final String allelesStr)
    {
        if(allelesStr == null || allelesStr.isEmpty())
            return Lists.newArrayList();

        String[] alleles = allelesStr.split(ITEM_DELIM, -1);
        return Arrays.stream(alleles).map(x -> HlaAllele.fromString(x)).collect(Collectors.toList());
    }
}
