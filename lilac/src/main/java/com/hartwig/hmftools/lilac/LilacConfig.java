package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import java.io.BufferedReader;
import java.io.InputStreamReader;
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

    public final String OutputFilePrefix;
    public final String OutputDir;
    public final String RefGenome;

    public final int MinBaseQual;
    public final int MinEvidence;
    public final int MinFragmentsPerAllele;
    public final int MinFragmentsToRemoveSingle;
    public final int MinConfirmedUniqueCoverage;
    public final int Threads;
    public final int MaxDistanceFromTopScore;

    public final String GeneCopyNumberFile;
    public final String SomaticVcf;

    public final boolean DebugPhasing;
    public final List<HlaAllele> ExpectedAlleles;
    public final List<HlaAllele> CommonAlleles;
    public final List<HlaAllele> StopLossRecoveryAlleles;

    // config strings
    private static final String SAMPLE = "sample";
    private static final String REFERENCE_BAM = "reference_bam";
    private static final String TUMOR_BAM = "tumor_bam";
    private static final String RESOURCE_DIR = "resource_dir";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String OUTPUT_ID = "output_ID";
    private static final String REF_GENOME = "ref_genome";
    private static final String MIN_BASE_QUAL = "min_base_qual";
    private static final String MIN_EVIDENCE = "min_evidence";
    private static final String THREADS = "threads";
    private static final String MIN_FRAGMENTS_PER_ALLELE = "min_fragments_per_allele";
    private static final String MIN_FRAGMENTS_TO_REMOVE_SINGLE = "min_fragments_to_remove_single";
    private static final String MIN_CONFIRMED_UNIQUE_COVERAGE = "min_confirmed_unique_coverage";
    private static final String EXPECTED_ALLELES = "expected_alleles";
    private static final String GENE_COPY_NUMBER = "gene_copy_number";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String MAX_DISTANCE_FROM_TOP_SCORE = "max_distance_from_top_score";

    private static final String DEBUG_PHASING = "debug_phasing";

    public static final Logger LL_LOGGER = LogManager.getLogger(LilacConfig.class);;

    public LilacConfig(final CommandLine cmd)
    {
        OutputDir = parseOutputDir(cmd);
        Sample = cmd.getOptionValue(SAMPLE);

        ReferenceBam = cmd.getOptionValue(REFERENCE_BAM);
        TumorBam = cmd.getOptionValue(TUMOR_BAM);
        ResourceDir = cmd.getOptionValue(RESOURCE_DIR);

        OutputFilePrefix = cmd.getOptionValue(OUTPUT_ID);
        RefGenome = cmd.getOptionValue(SAMPLE);

        MinBaseQual = ConfigUtils.getConfigValue(cmd, MIN_BASE_QUAL, LilacConstants.DEFAULT_MIN_BASE_QUAL);
        MinEvidence = ConfigUtils.getConfigValue(cmd, MIN_EVIDENCE, LilacConstants.DEFAULT_MIN_EVIDENCE);
        MinFragmentsPerAllele = ConfigUtils.getConfigValue(cmd, MIN_FRAGMENTS_PER_ALLELE, LilacConstants.DEFAULT_FRAGS_PER_ALLELE);
        MinFragmentsToRemoveSingle = ConfigUtils.getConfigValue(cmd, MIN_FRAGMENTS_TO_REMOVE_SINGLE, LilacConstants.DEFAULT_FRAGS_REMOVE_SGL);
        MinConfirmedUniqueCoverage = ConfigUtils.getConfigValue(cmd, MIN_CONFIRMED_UNIQUE_COVERAGE, LilacConstants.DEFAULT_MIN_CONF_UNIQUE_COVERAGE);

        MaxDistanceFromTopScore = ConfigUtils.getConfigValue(cmd, MAX_DISTANCE_FROM_TOP_SCORE, LilacConstants.DEFAULT_MAX_DIST_FROM_TOP_SCORE);

        GeneCopyNumberFile = cmd.getOptionValue(GENE_COPY_NUMBER);
        SomaticVcf = cmd.getOptionValue(SOMATIC_VCF);

        ExpectedAlleles = parseExpectedAlleles(cmd.getOptionValue(EXPECTED_ALLELES));
        CommonAlleles = loadCommonAlleles();
        StopLossRecoveryAlleles = loadStopLossRecoveryAllele();

        Threads = getConfigValue(cmd, THREADS, 1);
        DebugPhasing = cmd.hasOption(DEBUG_PHASING);
    }

    public LilacConfig()
    {
        OutputDir = "";
        Sample = "";
        ReferenceBam = "";
        TumorBam = "";
        ResourceDir = "";
        RefGenome = "";

        MinBaseQual = LilacConstants.DEFAULT_MIN_BASE_QUAL;
        MinEvidence = LilacConstants.DEFAULT_MIN_EVIDENCE;
        MinFragmentsPerAllele = LilacConstants.DEFAULT_FRAGS_PER_ALLELE;
        MinFragmentsToRemoveSingle = LilacConstants.DEFAULT_FRAGS_REMOVE_SGL;
        MinConfirmedUniqueCoverage = LilacConstants.DEFAULT_MIN_CONF_UNIQUE_COVERAGE;
        MaxDistanceFromTopScore = LilacConstants.DEFAULT_MAX_DIST_FROM_TOP_SCORE;

        GeneCopyNumberFile = "";
        SomaticVcf = "";

        // public final boolean DebugPhasing;
        ExpectedAlleles = Lists.newArrayList();
        CommonAlleles = loadCommonAlleles();
        StopLossRecoveryAlleles = loadStopLossRecoveryAllele();

        OutputFilePrefix = "";
        Threads = 0;
        DebugPhasing = false;
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(REFERENCE_BAM, true,"Path to reference/normal bam");
        options.addOption(TUMOR_BAM, true,"Path to tumor bam");
        options.addOption(RESOURCE_DIR, true,"Path to resource files");
        options.addOption(OUTPUT_DIR, true,"Path to output");
        options.addOption(OUTPUT_ID, true,"Output file prefix");
        options.addOption(REF_GENOME, true,"Optional path to reference genome fasta file");
        options.addOption(MIN_BASE_QUAL, true,"MIN_BASE_QUAL");
        options.addOption(MIN_EVIDENCE, true,"MIN_EVIDENCE");
        options.addOption(MIN_FRAGMENTS_PER_ALLELE, true,"MIN_FRAGMENTS_PER_ALLELE");
        options.addOption(MIN_FRAGMENTS_TO_REMOVE_SINGLE, true,"MIN_FRAGMENTS_TO_REMOVE_SINGLE");
        options.addOption(MIN_CONFIRMED_UNIQUE_COVERAGE, true,"MIN_CONFIRMED_UNIQUE_COVERAGE");
        options.addOption(MAX_DISTANCE_FROM_TOP_SCORE, true,"Max distance from top score");
        options.addOption(THREADS, true,"Number of threads");
        options.addOption(EXPECTED_ALLELES, true,"Comma separated expected alleles for the sample");
        options.addOption(GENE_COPY_NUMBER, true,"Path to gene copy number file");
        options.addOption(SOMATIC_VCF, true,"Path to somatic VCF");
        options.addOption(DEBUG_PHASING, false, "More detailed logging of phasing");
        DatabaseAccess.addDatabaseCmdLineArgs((Options) options);
        return options;
    }

    private final List<HlaAllele> parseExpectedAlleles(final String allelesStr)
    {
        if(allelesStr == null || allelesStr.isEmpty())
            return Lists.newArrayList();

        String[] alleles = allelesStr.split(";", -1);
        return Arrays.stream(alleles).map(x -> HlaAllele.fromString(x)).collect(Collectors.toList());
    }

    private final List<HlaAllele> loadCommonAlleles()
    {
        final List<String> commonAlleleLines = new BufferedReader(new InputStreamReader(
                RefGenomeCoordinates.class.getResourceAsStream("/alleles/common.txt")))
                .lines().collect(Collectors.toList());

        return commonAlleleLines.stream().map(x -> HlaAllele.fromString(x)).map(x -> x.asFourDigit()).collect(Collectors.toList());
    }

    private static List<HlaAllele> loadStopLossRecoveryAllele()
    {
        return Lists.newArrayList(HlaAllele.fromString("C*04:09N"));
    }

    /*
    @NotNull
    public final String requiredFile$lilac(@NotNull CommandLine $receiver, @NotNull String argument) throws IOException
    {
        Intrinsics.checkParameterIsNotNull((Object) $receiver, (String) "$receiver");
        Intrinsics.checkParameterIsNotNull((Object) argument, (String) "argument");
        String result = $receiver.getOptionValue(argument);
        if(!new File(result).exists())
        {
            throw (Throwable) new IOException("Unable to read file " + result);
        }
        String string = result;
        Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "result");
        return string;
    }

    @NotNull
    public final String optionalFile$lilac(@NotNull CommandLine $receiver, @NotNull String argument, @NotNull String string)
            throws IOException
    {
        Intrinsics.checkParameterIsNotNull((Object) $receiver, (String) "$receiver");
        Intrinsics.checkParameterIsNotNull((Object) argument, (String) "argument");
        Intrinsics.checkParameterIsNotNull((Object) string, (String) "default");
        if($receiver.hasOption(argument))
        {
            return this.requiredFile$lilac($receiver, argument);
        }
        return string;
    }
    */

    /*
    @NotNull
    public final String requiredDir$lilac(@NotNull CommandLine $receiver, @NotNull String argument) throws IOException
    {
        Intrinsics.checkParameterIsNotNull((Object) $receiver, (String) "$receiver");
        Intrinsics.checkParameterIsNotNull((Object) argument, (String) "argument");
        String result = $receiver.getOptionValue(argument);
        File dir = new File(result);
        if(!dir.exists() && !dir.mkdirs())
        {
            throw (Throwable) new IOException("Unable to create director " + result);
        }
        String string = result;
        Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "result");
        return string;
    }
     */

}
