package com.hartwig.hmftools.cobalt;

import static com.hartwig.hmftools.common.cli.Configs.defaultEnumValue;
import static com.hartwig.hmftools.common.cli.Configs.defaultIntValue;

import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cli.Configs;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.ValidationStringency;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltConfig {

    int DEFAULT_THREADS = 4;
    int DEFAULT_MIN_MAPPING_QUALITY = 10;

    String TUMOR_ONLY = "tumor_only";
    String TUMOR_ONLY_DIPLOID_BED = "tumor_only_diploid_bed";
    String THREADS = "threads";
    String REFERENCE = "reference";
    String REFERENCE_BAM = "reference_bam";
    String TUMOR = "tumor";
    String TUMOR_BAM = "tumor_bam";
    String REF_GENOME = "ref_genome";
    String OUTPUT_DIR = "output_dir";
    String INPUT_DIR = "input_dir";
    String GC_PROFILE = "gc_profile";
    String MIN_MAPPING_QUALITY = "min_quality";
    String VALIDATION_STRINGENCY = "validation_stringency";

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(TUMOR_ONLY, false, "Tumor only mode");
        options.addOption(TUMOR_ONLY_DIPLOID_BED, true, "Diploid regions for tumor-only mode");
        options.addOption(THREADS, true, "Number of threads [" + DEFAULT_THREADS + "]");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(INPUT_DIR, true, "Input directory (used for migration only)");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(MIN_MAPPING_QUALITY, true, "Min quality [" + DEFAULT_MIN_MAPPING_QUALITY + "]");
        options.addOption(GC_PROFILE, true, "Location of GC Profile");
        options.addOption(REF_GENOME, true, "Path to reference genome fasta file if using CRAM files");
        options.addOption(VALIDATION_STRINGENCY, true, "SAM validation strategy: STRICT, SILENT, LENIENT [STRICT]");

        return options;
    }

    int threadCount();

    int minMappingQuality();

    @NotNull
    String gcProfilePath();

    @NotNull
    String tumorBamPath();

    @NotNull
    String referenceBamPath();

    @NotNull
    String refGenomePath();

    @NotNull
    String inputDirectory();

    @NotNull
    String outputDirectory();

    @NotNull
    String reference();

    @NotNull
    String tumor();

    @NotNull
    ValidationStringency validationStringency();

    boolean tumorOnly();

    @NotNull
    String tumorOnlyDiploidBed();

    default int windowSize() {
        return 1000;
    }

    @NotNull
    static CobaltConfig createConfig(@NotNull final CommandLine cmd) throws ParseException, IOException {
        final int threadCount = defaultIntValue(cmd, THREADS, DEFAULT_THREADS);
        final int minMappingQuality = defaultIntValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY);
        final String refGenomePath = cmd.getOptionValue(REF_GENOME, "");

        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String gcProfilePath = Configs.readableFile(cmd, GC_PROFILE);
        if (gcProfilePath.endsWith("gz")) {
            throw new ParseException("Please supply un-compressed " + GC_PROFILE + " file");
        }

        if (cmd.hasOption(INPUT_DIR)) {
            throw new ParseException(INPUT_DIR + " not applicable to COBALT");
        }

        final boolean isTumorOnly = cmd.hasOption(TUMOR_ONLY);
        final String reference;
        final String referenceBamPath;
        final String diploidBed;
        if (isTumorOnly) {
            if (cmd.hasOption(REFERENCE_BAM)) {
                throw new ParseException(REFERENCE_BAM + " not applicable to tumor only mode");
            }

            if (cmd.hasOption(REFERENCE)) {
                throw new ParseException(REFERENCE + " not applicable to tumor only mode");
            }

            reference = CobaltRatioFile.TUMOR_ONLY_REFERENCE_SAMPLE;
            referenceBamPath = Strings.EMPTY;
            diploidBed = Configs.readableFile(cmd, TUMOR_ONLY_DIPLOID_BED);

        } else {
            diploidBed = Strings.EMPTY;
            reference = parameter(cmd, REFERENCE, missingJoiner);
            referenceBamPath = Configs.readableFile(cmd, REFERENCE_BAM);
        }

        final String tumorBamPath = Configs.readableFile(cmd, TUMOR_BAM);
        final String outputDirectory = Configs.writableOutputDirectory(cmd, OUTPUT_DIR);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String missing = missingJoiner.toString();

        final ValidationStringency validationStringency =
                defaultEnumValue(cmd, VALIDATION_STRINGENCY, ValidationStringency.DEFAULT_STRINGENCY);

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        return ImmutableCobaltConfig.builder()
                .threadCount(threadCount)
                .tumorOnly(isTumorOnly)
                .tumorOnlyDiploidBed(diploidBed)
                .minMappingQuality(minMappingQuality)
                .gcProfilePath(gcProfilePath)
                .tumorBamPath(tumorBamPath)
                .inputDirectory(Strings.EMPTY)
                .referenceBamPath(referenceBamPath)
                .refGenomePath(refGenomePath)
                .outputDirectory(outputDirectory)
                .reference(reference)
                .tumor(tumor)
                .validationStringency(validationStringency)
                .build();
    }

    @NotNull
    static CobaltConfig createMigrationConfig(@NotNull final CommandLine cmd) throws ParseException {
        final StringJoiner missingJoiner = new StringJoiner(", ");

        final String gcProfilePath = parameter(cmd, GC_PROFILE, missingJoiner);
        if (gcProfilePath.endsWith("gz")) {
            throw new ParseException("Please supply un-compressed " + GC_PROFILE + " file");
        }

        final String normal = parameter(cmd, REFERENCE, missingJoiner);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String inputDirectory = parameter(cmd, INPUT_DIR, missingJoiner);
        final String outputDirectory = parameter(cmd, OUTPUT_DIR, missingJoiner);
        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        if (inputDirectory.equals(outputDirectory)) {
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

    @NotNull
    static String parameter(@NotNull final CommandLine cmd, @NotNull final String parameter, @NotNull final StringJoiner missing) {
        final String value = cmd.getOptionValue(parameter);
        if (value == null) {
            missing.add(parameter);
            return "";
        }
        return value;
    }
}
