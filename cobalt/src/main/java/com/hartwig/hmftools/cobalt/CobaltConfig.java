package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_MIN_MAPPING_QUALITY;
import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.ValidationStringency;

public class CobaltConfig
{
    public final int ThreadCount;

    public final int MinMappingQuality;

    public final String ReferenceId;
    public final String TumorId;

    public final String TumorBamPath;
    public final String ReferenceBamPath;
    public final String RefGenomePath;
    public final String GcProfilePath;
    public final String OutputDir;

    public final ValidationStringency Stringency;
    public final boolean TumorOnly;
    public final String TumorOnlyDiploidBed;

    private static final int DEFAULT_THREADS = 4;

    public static String TUMOR = "tumor";
    public static String REFERENCE = "reference";

    private static String TUMOR_ONLY = "tumor_only";
    private static String TUMOR_ONLY_DIPLOID_BED = "tumor_only_diploid_bed";
    private static String THREADS = "threads";
    private static String REFERENCE_BAM = "reference_bam";
    private static String TUMOR_BAM = "tumor_bam";
    private static String REF_GENOME = "ref_genome";
    private static String GC_PROFILE = "gc_profile";
    private static String MIN_MAPPING_QUALITY = "min_quality";
    private static String VALIDATION_STRINGENCY = "validation_stringency";

    public static final Logger CB_LOGGER = LogManager.getLogger(CobaltConfig.class);

    public CobaltConfig(final CommandLine cmd)
    {
        ThreadCount = getConfigValue(cmd, THREADS, DEFAULT_THREADS);
        MinMappingQuality = getConfigValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY);
        RefGenomePath = cmd.getOptionValue(REF_GENOME, "");

        final StringJoiner missingJoiner = new StringJoiner(", ");
        GcProfilePath = cmd.getOptionValue(GC_PROFILE);

        TumorOnly = cmd.hasOption(TUMOR_ONLY);

        if(TumorOnly)
        {
            ReferenceId = CobaltRatioFile.TUMOR_ONLY_REFERENCE_SAMPLE;
            ReferenceBamPath = Strings.EMPTY;
            TumorOnlyDiploidBed = cmd.getOptionValue(TUMOR_ONLY_DIPLOID_BED);
        }
        else
        {
            TumorOnlyDiploidBed = "";
            ReferenceId = parameter(cmd, REFERENCE, missingJoiner);
            ReferenceBamPath = cmd.getOptionValue(REFERENCE_BAM);
        }

        TumorBamPath = cmd.getOptionValue(TUMOR_BAM);
        OutputDir = parseOutputDir(cmd);

        TumorId = parameter(cmd, TUMOR, missingJoiner);

        Stringency = ValidationStringency.valueOf(cmd.getOptionValue(VALIDATION_STRINGENCY, "DEFAULT_STRINGENCY"));
    }

    public boolean isValid()
    {
        if(TumorOnly)
        {
            if(ReferenceBamPath != null || ReferenceId != null)
            {
                CB_LOGGER.error( "referenceId or BAM not applicable to tumor only mode");
                return false;
            }
        }

        if(!checkCreateOutputDir(OutputDir))
        {
            CB_LOGGER.error("failed to create output directory({})", OutputDir);
            return false;
        }

        if(GcProfilePath.endsWith("gz"))
        {
            CB_LOGGER.error( "invalid GC-profile file, must be uncompressed");
            return false;
        }

        return true;
    }

    public static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(TUMOR_ONLY, false, "Tumor only mode");
        options.addOption(TUMOR_ONLY_DIPLOID_BED, true, "Diploid regions for tumor-only mode");
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(MIN_MAPPING_QUALITY, true, "Min quality [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(GC_PROFILE, true, "Location of GC Profile");
        options.addOption(REF_GENOME, true, "Path to reference genome fasta file if using CRAM files");
        options.addOption(VALIDATION_STRINGENCY, true, "SAM validation strategy: STRICT, SILENT, LENIENT [STRICT]");

        return options;
    }

    /*
    public CobaltConfig createMigrationConfig(final CommandLine cmd) throws ParseException
    {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        GcProfilePath = parameter(cmd, GC_PROFILE, missingJoiner);
        if(gcProfilePath.endsWith("gz"))
        {
            throw new ParseException("Please supply un-compressed " + GC_PROFILE + " file");
        }

        final String normal = parameter(cmd, REFERENCE, missingJoiner);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String inputDirectory = parameter(cmd, INPUT_DIR, missingJoiner);
        final String outputDirectory = parameter(cmd, OUTPUT_DIR, missingJoiner);
        final String missing = missingJoiner.toString();

        if(!missing.isEmpty())
        {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        if(inputDirectory.equals(outputDirectory))
        {
            throw new ParseException("Input and output directories must be difference");
        }

        return ImmutableCobaltConfig.builder()
                .threadCount(1)
                .minMappingQuality(0)
                .gcProfilePath(gcProfilePath)
                .tumorBamPath(Strings.EMPTY)
                .referenceBamPath(Strings.EMPTY)
                .refGenomePath(Strings.EMPTY)
                .inputDirectory(inputDirectory)
                .outputDirectory(outputDirectory)
                .reference(normal)
                .tumor(tumor)
                .validationStringency(ValidationStringency.DEFAULT_STRINGENCY)
                .build();
    }
    */

    private String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing)
    {
        final String value = cmd.getOptionValue(parameter);
        if(value == null)
        {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
