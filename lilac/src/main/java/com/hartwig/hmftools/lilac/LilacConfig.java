package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FRAGS_PER_ALLELE;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_FRAGS_REMOVE_SGL;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_BASE_QUAL;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_MIN_EVIDENCE;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFAULT_TOP_SCORE_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.ITEM_DELIM;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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
    public final String Sample;
    public final String ReferenceBam;
    public final String TumorBam;
    public final String ResourceDir;

    public final String OutputDir;
    public final String RefGenome;

    public final int MinBaseQual;
    public final int MinEvidence;
    public final int MinFragmentsPerAllele;
    public final int MinFragmentsToRemoveSingle;
    public final int Threads;
    public final double TopScoreThreshold;

    public final String GeneCopyNumberFile;
    public final String SomaticVcf;

    public final boolean DebugPhasing;
    public final boolean RunValidation;

    // pre-determined sample alleles, always kept and scored, primarily for testing
    public final List<HlaAllele> ExpectedAlleles;

    // limit analysis from full data set to these + other configured alleles, for testing purposes
    public final List<HlaAllele> RestrictedAlleles;

    // config strings
    public static final String RESOURCE_DIR = "resource_dir";
    public static final String OUTPUT_DIR = "output_dir";

    private static final String SAMPLE = "sample";
    private static final String SAMPLE_DATA_DIR = "sample_data_dir";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String REF_GENOME = "ref_genome";
    private static final String MIN_BASE_QUAL = "min_base_qual";
    private static final String MIN_EVIDENCE = "min_evidence";
    private static final String THREADS = "threads";
    private static final String MIN_FRAGMENTS_PER_ALLELE = "min_fragments_per_allele";
    private static final String MIN_FRAGMENTS_TO_REMOVE_SINGLE = "min_fragments_to_remove_single";
    private static final String EXPECTED_ALLELES = "expected_alleles";
    private static final String RESTRICTED_ALLELES = "restricted_alleles";
    private static final String GENE_COPY_NUMBER = "gene_copy_number";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String TOP_SCORE_THRESHOLD = "top_score_threshold";

    private static final String DEBUG_PHASING = "debug_phasing";
    public static final String LOG_DEBUG = "log_debug";
    public static final String RUN_VALIDATION = "run_validation";
    public static final String LOG_LEVEL = "log_level";

    public static final Logger LL_LOGGER = LogManager.getLogger(LilacConfig.class);;

    public LilacConfig(final CommandLine cmd)
    {
        Sample = cmd.getOptionValue(SAMPLE);

        if(cmd.hasOption(SAMPLE_DATA_DIR))
        {
            String sampleDataDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DATA_DIR));

            OutputDir = sampleDataDir;

            String referenceId = Sample.substring(0, Sample.lastIndexOf('T')) + "R";
            ReferenceBam = sampleDataDir + referenceId + ".hla.bam";

            // these 3 are optional so check existence before initialising
            TumorBam = checkFileExists(sampleDataDir + Sample + ".hla.bam");
            GeneCopyNumberFile = checkFileExists(sampleDataDir + Sample + ".purple.cnv.gene.tsv");
            SomaticVcf = checkFileExists(sampleDataDir + Sample + ".purple.somatic.vcf.gz");
        }
        else
        {
            ReferenceBam = cmd.getOptionValue(REFERENCE_BAM, "");
            TumorBam = cmd.getOptionValue(TUMOR_BAM, "");
            GeneCopyNumberFile = cmd.getOptionValue(GENE_COPY_NUMBER, "");
            SomaticVcf = cmd.getOptionValue(SOMATIC_VCF, "");
            OutputDir = parseOutputDir(cmd);
        }

        ResourceDir = checkAddDirSeparator(cmd.getOptionValue(RESOURCE_DIR));
        RefGenome = cmd.getOptionValue(REF_GENOME, "");

        MinBaseQual = ConfigUtils.getConfigValue(cmd, MIN_BASE_QUAL, DEFAULT_MIN_BASE_QUAL);
        MinEvidence = ConfigUtils.getConfigValue(cmd, MIN_EVIDENCE, DEFAULT_MIN_EVIDENCE);
        MinFragmentsPerAllele = ConfigUtils.getConfigValue(cmd, MIN_FRAGMENTS_PER_ALLELE, DEFAULT_FRAGS_PER_ALLELE);
        MinFragmentsToRemoveSingle = ConfigUtils.getConfigValue(cmd, MIN_FRAGMENTS_TO_REMOVE_SINGLE, DEFAULT_FRAGS_REMOVE_SGL);

        TopScoreThreshold = ConfigUtils.getConfigValue(cmd, TOP_SCORE_THRESHOLD, DEFAULT_TOP_SCORE_THRESHOLD);

        ExpectedAlleles = parseAlleleList(cmd.getOptionValue(EXPECTED_ALLELES));
        RestrictedAlleles = parseAlleleList(cmd.getOptionValue(RESTRICTED_ALLELES));

        Threads = getConfigValue(cmd, THREADS, 1);

        DebugPhasing = cmd.hasOption(DEBUG_PHASING); //  || cmd.hasOption(LOG_);
        RunValidation = cmd.hasOption(RUN_VALIDATION);
    }

    private String checkFileExists(final String filename)
    {
        return Files.exists(Paths.get(filename)) ? filename : "";
    }

    public String outputPrefix() { return OutputDir + Sample;}

    public boolean isValid()
    {
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
        LL_LOGGER.info("sample({}) inputs: tumorBam({}) somaticVCF({}) geneCopyNumber({})",
                Sample, !TumorBam.isEmpty(), !SomaticVcf.isEmpty(), !GeneCopyNumberFile.isEmpty());

        LL_LOGGER.info("minBaseQual({}), minEvidence({}) minFragmentsPerAllele({}) "
                + "minFragmentsToRemoveSingle({}) maxDistanceFromTopScore({})",
                MinBaseQual, MinEvidence, MinFragmentsPerAllele, MinFragmentsToRemoveSingle, TopScoreThreshold);

        if(RunValidation)
        {
            LL_LOGGER.info("running validation routines");
        }

        if(!ExpectedAlleles.isEmpty())
        {
            LL_LOGGER.info("expected alleles: {}", HlaAllele.toString(ExpectedAlleles));
        }

        if(!RestrictedAlleles.isEmpty())
        {
            LL_LOGGER.info("restricted alleles: {}", HlaAllele.toString(RestrictedAlleles));
        }
    }

    public LilacConfig()
    {
        OutputDir = "";
        Sample = "";
        ReferenceBam = "";
        TumorBam = "";
        ResourceDir = "";
        RefGenome = "";

        MinBaseQual = DEFAULT_MIN_BASE_QUAL;
        MinEvidence = DEFAULT_MIN_EVIDENCE;
        MinFragmentsPerAllele = DEFAULT_FRAGS_PER_ALLELE;
        MinFragmentsToRemoveSingle = DEFAULT_FRAGS_REMOVE_SGL;
        TopScoreThreshold = DEFAULT_TOP_SCORE_THRESHOLD;

        GeneCopyNumberFile = "";
        SomaticVcf = "";

        ExpectedAlleles = Lists.newArrayList();
        RestrictedAlleles = Lists.newArrayList();

        Threads = 0;
        DebugPhasing = false;
        RunValidation = true;
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(SAMPLE_DATA_DIR, true,"Path to all sample files");
        options.addOption(REFERENCE_BAM, true,"Path to reference/normal bam");
        options.addOption(TUMOR_BAM, true,"Path to tumor bam");
        options.addOption(RESOURCE_DIR, true,"Path to resource files");
        options.addOption(OUTPUT_DIR, true,"Path to output");
        options.addOption(REF_GENOME, true,"Optional path to reference genome fasta file");
        options.addOption(MIN_BASE_QUAL, true,"MIN_BASE_QUAL");
        options.addOption(MIN_EVIDENCE, true,"MIN_EVIDENCE");
        options.addOption(MIN_FRAGMENTS_PER_ALLELE, true,"MIN_FRAGMENTS_PER_ALLELE");
        options.addOption(MIN_FRAGMENTS_TO_REMOVE_SINGLE, true,"MIN_FRAGMENTS_TO_REMOVE_SINGLE");
        options.addOption(TOP_SCORE_THRESHOLD, true,"Max distance from top score");
        options.addOption(THREADS, true,"Number of threads");
        options.addOption(EXPECTED_ALLELES, true,"Comma separated expected alleles for the sample");
        options.addOption(RESTRICTED_ALLELES, true,"Comma separated restricted analysis allele list");
        options.addOption(GENE_COPY_NUMBER, true,"Path to gene copy number file");
        options.addOption(SOMATIC_VCF, true,"Path to somatic VCF");
        options.addOption(DEBUG_PHASING, false, "More detailed logging of phasing");
        options.addOption(RUN_VALIDATION, false, "Run validation checks");
        options.addOption(LOG_DEBUG, false, "More detailed logging of phasing");
        options.addOption(LOG_LEVEL, true, "Specify log level WARN, INFO, DEBUG or TRACE");
        DatabaseAccess.addDatabaseCmdLineArgs((Options) options);
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
