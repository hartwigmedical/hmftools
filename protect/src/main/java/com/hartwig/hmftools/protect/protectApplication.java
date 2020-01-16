package com.hartwig.hmftools.protect;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class protectApplication {
    private static final Logger LOGGER = LogManager.getLogger(protectApplication.class);

    private static final String TUMOR_SAMPLE_ID = "tumor_sample_id";

    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";

    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        String tumorSampleId = cmd.getOptionValue(TUMOR_SAMPLE_ID);

        // General params needed for every sample
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);

        // Params specific for specific sample
        final String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        final String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);
        final String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        final String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);

        if (!validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        LOGGER.info("Create actionability for sample {}", tumorSampleId);

        LOGGER.info("Create actionability for patient report");

        LOGGER.info("Create actionability for database");

    }

    private static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, TUMOR_SAMPLE_ID) && dirExists(cmd, KNOWLEDGEBASE_DIRECTORY) && fileExists(cmd, TUMOR_LOCATION_CSV) && fileExists(cmd,
                SOMATIC_VARIANT_VCF) && fileExists(cmd, PURPLE_PURITY_TSV)
                && fileExists(cmd, PURPLE_GENE_CNV_TSV) && fileExists(cmd, LINX_FUSION_TSV);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            LOGGER.warn(param + " has to be provided");
            return false;
        }
        return true;
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn(param + " has to be an existing file: " + value);
            return false;
        }

        return true;
    }

    private static boolean dirExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value) || !pathIsDirectory(value)) {
            LOGGER.warn(param + " has to be an existing directory: " + value);
            return false;
        }

        return true;
    }

    private static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a patient report will be generated.");

        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the folder containing knowledgebase files.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");

        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");

        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Protect", options);
        System.exit(1);
    }

}
