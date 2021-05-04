package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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
    public final List<com.hartwig.hmftools.lilac.hla.HlaAllele> ExpectedAlleles;
    public final List<com.hartwig.hmftools.lilac.hla.HlaAllele> CommonAlleles;
    public final List<HlaAllele> StopLossRecoveryAlleles;

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

    public LilacConfig(final CommandLine cmd) throws ParseException, IOException
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

        // public final boolean DebugPhasing;
        ExpectedAlleles = Lists.newArrayList();
        CommonAlleles = Lists.newArrayList();
        StopLossRecoveryAlleles = Lists.newArrayList();

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
        CommonAlleles = Lists.newArrayList();
        StopLossRecoveryAlleles = Lists.newArrayList();

        /*
        val expectedAlleles = cmd.expectedAlleles(EXPECTED_ALLELES);
        val commonAlleles = LilacConfig::class.java.getResource("/alleles/common.txt")
            .readText()
            .split("\n")
            .map { HlaAllele(it) }
                    .map { it.asFourDigit() }

        val stopLossRecoveryAlleles = listOf(HlaAllele("C*04:09N"))
        */

        OutputFilePrefix = "";
        Threads = 0;
        DebugPhasing = false;
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, "Name of sample");
        options.addOption(REFERENCE_BAM, "Path to reference/normal bam");
        options.addOption(TUMOR_BAM, "Path to tumor bam");
        options.addOption(RESOURCE_DIR, "Path to resource files");
        options.addOption(OUTPUT_DIR, "Path to output");
        options.addOption(OUTPUT_ID, "Output file prefix");
        options.addOption(REF_GENOME, "Optional path to reference genome fasta file");
        options.addOption(MIN_BASE_QUAL, "MIN_BASE_QUAL");
        options.addOption(MIN_EVIDENCE, "MIN_EVIDENCE");
        options.addOption(MIN_FRAGMENTS_PER_ALLELE, "MIN_FRAGMENTS_PER_ALLELE");
        options.addOption(MIN_FRAGMENTS_TO_REMOVE_SINGLE, "MIN_FRAGMENTS_TO_REMOVE_SINGLE");
        options.addOption(MIN_CONFIRMED_UNIQUE_COVERAGE, "MIN_CONFIRMED_UNIQUE_COVERAGE");
        options.addOption(MAX_DISTANCE_FROM_TOP_SCORE, "Max distance from top score");
        options.addOption(THREADS, "Number of threads");
        options.addOption(EXPECTED_ALLELES, "Comma separated expected alleles");
        options.addOption(GENE_COPY_NUMBER, "Path to gene copy number file");
        options.addOption(SOMATIC_VCF, "Path to somatic VCF");
        options.addOption(DEBUG_PHASING, false, "More detailed logging of phasing");
        DatabaseAccess.addDatabaseCmdLineArgs((Options) options);
        return options;
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
    private final List<HlaAllele> expectedAlleles(CommandLine $receiver, String opt)
    {
        if($receiver.hasOption(opt))
        {
            void $receiver$iv$iv;
            Iterable $receiver$iv;
            String string = $receiver.getOptionValue(opt);
            if(string == null)
            {
                Intrinsics.throwNpe();
            }
            Iterable iterable = $receiver$iv = (Iterable) StringsKt.split$default((CharSequence) string, (String[]) new String[] {
                    "," }, (boolean) false, (int) 0, (int) 6, null);
            Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
            for(Object item$iv$iv : $receiver$iv$iv)
            {
                void it;
                String string2 = (String) item$iv$iv;
                Collection collection = destination$iv$iv;
                boolean bl = false;
                HlaAllele hlaAllele = HlaAllele.Companion.invoke((String) it).asFourDigit();
                collection.add(hlaAllele);
            }
            List result = (List) destination$iv$iv;
            this.getLogger()
                    .info("Using non default value {} for parameter {}", (Object) CollectionsKt.joinToString$default((Iterable) result, (CharSequence) ",", null, null, (int) 0, null, null, (int) 62, null), (Object) opt);
            return result;
        }
        return CollectionsKt.emptyList();
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
