package com.hartwig.hmftools.sage;

import java.util.StringJoiner;

import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface SageHotspotApplicationConfig {

    Logger LOGGER = LogManager.getLogger(SageHotspotApplicationConfig.class);

    String OUT_PATH = "out";
    String TUMOR = "tumor";
    String TUMOR_BAM = "tumor_bam";
    String REFERENCE = "reference";
    String REF_GENOME = "ref_genome";
    String REFERENCE_BAM = "reference_bam";
    String CODING_REGIONS = "coding_regions";
    String KNOWN_HOTSPOTS = "known_hotspots";

    String MIN_SNV_VAF = "min_snv_vaf";
    String MIN_INDEL_VAF = "min_indel_vaf";
    String MIN_SNV_QUALITY = "min_snv_quality";
    String MIN_INDEL_QUALITY = "min_indel_quality";

    String MIN_TUMOR_READS = "min_tumor_reads";
    String MIN_BASE_QUALITY = "min_base_quality";
    String MIN_MAPPING_QUALITY = "min_mapping_quality";
    String MAX_HET_BINOMIAL_LIKELIHOOD = "max_het_binomial_likelihood";

    int DEFAULT_MIN_TUMOR_READS = 2;
    int DEFAULT_MIN_BASE_QUALITY = 13;
    int DEFAULT_MIN_MAPPING_QUALITY = 1;
    int DEFAULT_MIN_SNV_QUALITY = 100;
    int DEFAULT_MIN_INDEL_QUALITY = 150;
    int DEFAULT_TYPICAL_READ_LENGTH = 151;

    double DEFAULT_MIN_SNV_VAF = 0.005;
    double DEFAULT_MIN_INDEL_VAF = 0.02;
    double DEFAULT_MAX_HET_BINOMIAL_LIKELIHOOD = 0.01;

    @NotNull
    static Options createOptions() {
        final Options options = new Options();
        options.addOption(OUT_PATH, true, "Output VCF file to write. Gz supported.");
        options.addOption(TUMOR, true, "Name of tumor sample");
        options.addOption(TUMOR_BAM, true, "Path to tumor bam file");
        options.addOption(REFERENCE, true, "Name of reference sample");
        options.addOption(REFERENCE_BAM, true, "Path to reference bam file");
        options.addOption(KNOWN_HOTSPOTS, true, "Tab separated file of known hotspot locations");
        options.addOption(CODING_REGIONS, true, "Coding regions bed file to search for inframe indels");
        options.addOption(MIN_BASE_QUALITY, true, "Minimum quality for a base to be considered [" + DEFAULT_MIN_BASE_QUALITY + "]");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");

        options.addOption(MAX_HET_BINOMIAL_LIKELIHOOD,
                true,
                "Up to 1 read of support for the ALT in the reference is allowed if the binomial likelihood of observing the read (assuming a 50% VAF) is less than this parameter ["
                        + DEFAULT_MAX_HET_BINOMIAL_LIKELIHOOD + "]");

        options.addOption(MIN_MAPPING_QUALITY,
                true,
                "Minimum mapping quality for an alignment to be used [" + DEFAULT_MIN_MAPPING_QUALITY + "]");

        options.addOption(MIN_TUMOR_READS, true, "Low confidence filtering minimum tumor reads [" + DEFAULT_MIN_TUMOR_READS + "]");
        options.addOption(MIN_SNV_VAF,
                true,
                "Low confidence filtering minimum SNV/MNV VAF [" + DEFAULT_MIN_SNV_VAF + "]");
        options.addOption(MIN_INDEL_VAF,
                true,
                "Low confidence filtering minimum indel VAF [" + DEFAULT_MIN_INDEL_VAF + "]");
        options.addOption(MIN_SNV_QUALITY,
                true,
                "Low confidence filtering minimum SNV/MNV quality [" + DEFAULT_MIN_SNV_QUALITY + "]");
        options.addOption(MIN_INDEL_QUALITY,
                true,
                "Low confidence filtering minimum indel quality [" + DEFAULT_MIN_INDEL_QUALITY + "]");

        return options;
    }

    @NotNull
    String codingRegionBedPath();

    @NotNull
    String tumorBamPath();

    @NotNull
    String referenceBamPath();

    @NotNull
    String refGenomePath();

    @NotNull
    String outputFile();

    @NotNull
    String normal();

    @NotNull
    String tumor();

    @NotNull
    String knownHotspotPath();

    int minTumorReads();

    int minSnvQuality();

    int minIndelQuality();

    int minMappingQuality();

    int minBaseQuality();

    double maxHetBinomialLikelihood();

    double minSnvVAF();

    double minIndelVAF();

    @NotNull
    static SageHotspotApplicationConfig createConfig(@NotNull final CommandLine cmd) throws ParseException {
        final StringJoiner missingJoiner = new StringJoiner(", ");
        final String codingRegions = parameter(cmd, CODING_REGIONS, missingJoiner);
        final String tumorBamPath = parameter(cmd, TUMOR_BAM, missingJoiner);
        final String referenceBamPath = parameter(cmd, REFERENCE_BAM, missingJoiner);
        final String refGenomePath = parameter(cmd, REF_GENOME, missingJoiner);
        final String outputVCF = parameter(cmd, OUT_PATH, missingJoiner);
        final String normal = parameter(cmd, REFERENCE, missingJoiner);
        final String tumor = parameter(cmd, TUMOR, missingJoiner);
        final String knownHotspotPath = parameter(cmd, KNOWN_HOTSPOTS, missingJoiner);
        final String missing = missingJoiner.toString();

        if (!missing.isEmpty()) {
            throw new ParseException("Missing the following parameters: " + missing);
        }

        return ImmutableSageHotspotApplicationConfig.builder()
                .tumor(tumor)
                .normal(normal)
                .outputFile(outputVCF)
                .tumorBamPath(tumorBamPath)
                .refGenomePath(refGenomePath)
                .codingRegionBedPath(codingRegions)
                .referenceBamPath(referenceBamPath)
                .knownHotspotPath(knownHotspotPath)
                .minTumorReads(defaultIntValue(cmd, MIN_TUMOR_READS, DEFAULT_MIN_TUMOR_READS))
                .minBaseQuality(defaultIntValue(cmd, MIN_BASE_QUALITY, DEFAULT_MIN_BASE_QUALITY))
                .minSnvVAF(defaultDoubleValue(cmd, MIN_SNV_VAF, DEFAULT_MIN_SNV_VAF))
                .minIndelVAF(defaultDoubleValue(cmd, MIN_INDEL_VAF, DEFAULT_MIN_INDEL_VAF))
                .minMappingQuality(defaultIntValue(cmd, MIN_MAPPING_QUALITY, DEFAULT_MIN_MAPPING_QUALITY))
                .minSnvQuality(defaultIntValue(cmd, MIN_SNV_QUALITY, DEFAULT_MIN_SNV_QUALITY))
                .minIndelQuality(defaultIntValue(cmd, MIN_INDEL_QUALITY, DEFAULT_MIN_INDEL_QUALITY))
                .maxHetBinomialLikelihood(defaultDoubleValue(cmd, MAX_HET_BINOMIAL_LIKELIHOOD, DEFAULT_MAX_HET_BINOMIAL_LIKELIHOOD))
                .build();
    }

    static double defaultDoubleValue(@NotNull final CommandLine cmd, @NotNull final String opt, final double defaultValue) {
        if (cmd.hasOption(opt)) {
            final double result = Double.parseDouble(cmd.getOptionValue(opt));
            if (!Doubles.equal(result, defaultValue)) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
    }

    static int defaultIntValue(@NotNull final CommandLine cmd, @NotNull final String opt, final int defaultValue) {
        if (cmd.hasOption(opt)) {
            final int result = Integer.parseInt(cmd.getOptionValue(opt));
            if (result != defaultValue) {
                LOGGER.info("Using non default value {} for parameter {}", result, opt);
            }
            return result;
        }

        return defaultValue;
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
