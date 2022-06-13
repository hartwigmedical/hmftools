package com.hartwig.hmftools.rose;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface RoseConfig {

    String OUTPUT_DIRECTORY = "output_dir";
    String ACTIONABILITY_DATABASE_TSV = "actionability_database_tsv";
    String TUMOR_SAMPLE_ID = "tumor_sample_id";
    String REF_SAMPLE_ID = "ref_sample_id";
    String PATIENT_ID = "patient_id";
    String PURPLE_PURITY_TSV = "purple_purity_tsv";
    String PURPLE_QC_FILE = "purple_qc_file";
    String PURPLE_SOMATIC_DRIVER_CATALOG_TSV = "purple_somatic_driver_catalog_tsv";
    String PURPLE_GERMLINE_DRIVER_CATALOG_TSV = "purple_germline_driver_catalog_tsv";
    String PURPLE_SOMATIC_VARIANT_VCF = "purple_somatic_variant_vcf";
    String PURPLE_GERMLINE_VARIANT_VCF = "purple_germline_variant_vcf";
    String PURPLE_GENE_COPY_NUMBER_TSV = "purple_gene_copy_number_tsv";
    String LINX_FUSION_TSV = "linx_fusion_tsv";
    String LINX_BREAKEND_TSV = "linx_breakend_tsv";
    String LINX_DRIVER_CATALOG_TSV = "linx_driver_catalog_tsv";
    String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    String ANNOTATED_VIRUS_TSV = "annotated_virus_tsv";
    String DRIVER_GENE_37_TSV = "driver_gene_37_tsv";
    String DRIVER_GENE_38_TSV = "driver_gene_38_tsv";
    String PRIMARY_TUMOR_TSV = "primary_tumor_tsv";
    String MOLECULAR_TISSUE_ORIGIN_TXT = "molecular_tissue_origin_txt";
    // Some additional optional params and flags
    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(OUTPUT_DIRECTORY, true, "Path to where the data of the report will be written to.");
        options.addOption(ACTIONABILITY_DATABASE_TSV, true, "Path to where the data oof the actionability database can be found.");
        options.addOption(RefGenomeVersion.REF_GENOME_VERSION, true, "Ref genome version to use (either '37' or '38')");

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a conclusion will be generated.");
        options.addOption(REF_SAMPLE_ID, true, "The reference ID for which is used for this sample.");
        options.addOption(PATIENT_ID, true, "The patient ID of the sample ID.");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");
        options.addOption(PURPLE_SOMATIC_DRIVER_CATALOG_TSV, true, "Path towards the purple somatic driver catalog TSV.");
        options.addOption(PURPLE_GERMLINE_DRIVER_CATALOG_TSV, true, "Path towards the purple germline driver catalog TSV.");
        options.addOption(PURPLE_SOMATIC_VARIANT_VCF, true, "Path towards the purple somatic variant VCF.");
        options.addOption(PURPLE_GERMLINE_VARIANT_VCF, true, "Path towards the purple germline variant VCF.");
        options.addOption(PURPLE_GENE_COPY_NUMBER_TSV, true, "Path towards the purple somatic copynumber TSV.");

        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_BREAKEND_TSV, true, "Path towards the linx breakend TSV.");
        options.addOption(LINX_DRIVER_CATALOG_TSV, true, "Path towards the LINX driver catalog TSV.");

        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT.");

        options.addOption(ANNOTATED_VIRUS_TSV, true, "Path towards the annotated virus TSV.");

        options.addOption(MOLECULAR_TISSUE_ORIGIN_TXT, true, "Path towards the molecular tissue origin TXT.");

        options.addOption(DRIVER_GENE_37_TSV, true, "Path to driver gene v37 TSV");
        options.addOption(DRIVER_GENE_38_TSV, true, "Path to driver gene v38 TSV");

        options.addOption(PRIMARY_TUMOR_TSV, true, "Path towards the (curated) primary tumor TSV.");

        options.addOption(LOG_DEBUG, false, "If provided, set the log level to debug rather than default.");
        return options;
    }

    @NotNull
    String outputDir();

    @NotNull
    String actionabilityDatabaseTsv();

    @NotNull
    RefGenomeVersion refGenomeVersion();

    @NotNull
    String tumorSampleId();

    @Nullable
    String refSampleId();

    @NotNull
    String patientId();

    @NotNull
    String purplePurityTsv();

    @NotNull
    String purpleQcFile();

    @NotNull
    String purpleSomaticDriverCatalogTsv();

    @NotNull
    String purpleGermlineDriverCatalogTsv();

    @NotNull
    String purpleSomaticVariantVcf();

    @NotNull
    String purpleGermlineVariantVcf();

    @NotNull
    String purpleSomaticCopyNumberTsv();

    @NotNull
    String linxFusionTsv();

    @NotNull
    String linxBreakendTsv();

    @NotNull
    String linxDriverCatalogTsv();

    @NotNull
    String chordPredictionTxt();

    @NotNull
    String annotatedVirusTsv();

    @NotNull
    String molecularTissueOriginTxt();

    @NotNull
    String driverGene37Tsv();

    @NotNull
    String driverGene38Tsv();

    @NotNull
    String primaryTumorTsv();

    @NotNull
    static RoseConfig createConfig(@NotNull CommandLine cmd) throws ParseException, IOException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        return ImmutableRoseConfig.builder()
                .outputDir(nonOptionalDir(cmd, OUTPUT_DIRECTORY))
                .actionabilityDatabaseTsv(nonOptionalFile(cmd, ACTIONABILITY_DATABASE_TSV))
                .refGenomeVersion(RefGenomeVersion.from(nonOptionalValue(cmd, RefGenomeVersion.REF_GENOME_VERSION)))
                .tumorSampleId(nonOptionalValue(cmd, TUMOR_SAMPLE_ID))
                .refSampleId(cmd.hasOption(REF_SAMPLE_ID) ? nonOptionalValue(cmd, REF_SAMPLE_ID) : null)
                .patientId(nonOptionalValue(cmd, PATIENT_ID))
                .purplePurityTsv(nonOptionalFile(cmd, PURPLE_PURITY_TSV))
                .purpleQcFile(nonOptionalFile(cmd, PURPLE_QC_FILE))
                .purpleSomaticDriverCatalogTsv(nonOptionalFile(cmd, PURPLE_SOMATIC_DRIVER_CATALOG_TSV))
                .purpleGermlineDriverCatalogTsv(nonOptionalFile(cmd, PURPLE_GERMLINE_DRIVER_CATALOG_TSV))
                .purpleSomaticVariantVcf(nonOptionalFile(cmd, PURPLE_SOMATIC_VARIANT_VCF))
                .purpleGermlineVariantVcf(nonOptionalFile(cmd, PURPLE_GERMLINE_VARIANT_VCF))
                .purpleSomaticCopyNumberTsv(nonOptionalFile(cmd, PURPLE_GENE_COPY_NUMBER_TSV))
                .linxFusionTsv(nonOptionalFile(cmd, LINX_FUSION_TSV))
                .linxBreakendTsv(nonOptionalFile(cmd, LINX_BREAKEND_TSV))
                .linxDriverCatalogTsv(nonOptionalFile(cmd, LINX_DRIVER_CATALOG_TSV))
                .chordPredictionTxt(nonOptionalFile(cmd, CHORD_PREDICTION_TXT))
                .annotatedVirusTsv(nonOptionalFile(cmd, ANNOTATED_VIRUS_TSV))
                .molecularTissueOriginTxt(nonOptionalFile(cmd, MOLECULAR_TISSUE_ORIGIN_TXT))
                .driverGene37Tsv(nonOptionalFile(cmd, DRIVER_GENE_37_TSV))
                .driverGene38Tsv(nonOptionalFile(cmd, DRIVER_GENE_38_TSV))
                .primaryTumorTsv(nonOptionalFile(cmd, PRIMARY_TUMOR_TSV))
                .build();
    }

    @NotNull
    static String nonOptionalValue(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            throw new ParseException("Parameter must be provided: " + param);
        }

        return value;
    }

    @NotNull
    static String nonOptionalDir(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value) || !pathIsDirectory(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing directory: " + value);
        }

        return value;
    }

    @NotNull
    static String nonOptionalFile(@NotNull CommandLine cmd, @NotNull String param) throws ParseException {
        String value = nonOptionalValue(cmd, param);

        if (!pathExists(value)) {
            throw new ParseException("Parameter '" + param + "' must be an existing file: " + value);
        }

        return value;
    }

    static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    static boolean pathIsDirectory(@NotNull String path) {
        return Files.isDirectory(new File(path).toPath());
    }
}