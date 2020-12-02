package com.hartwig.hmftools.serve;

import java.io.File;
import java.nio.file.Files;

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
public interface ServeConfig {

    // Input sources to SERVE
    String VICC_JSON = "vicc_json";
    String ICLUSION_TRIAL_TSV = "iclusion_trial_tsv";
    String DOCM_TSV = "docm_tsv";
    String HARTWIG_COHORT_TSV = "hartwig_cohort_tsv";
    String HARTWIG_CURATED_TSV = "hartwig_curated_tsv";

    // Additional config for knowledge generation
    String REF_GENOME_VERSION = "ref_genome_version";
    String REF_GENOME_FASTA_FILE = "ref_genome_fasta_file";
    String DRIVER_GENE_TSV = "driver_gene_tsv";

    // All output from SERVE will be written to this dir
    String OUTPUT_DIR = "output_dir";

    // Options to help with debugging / testing
    String SKIP_HOTSPOT_RESOLVING = "skip_hotspot_resolving";
    String LOG_DEBUG = "log_debug";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(VICC_JSON, true, "Path to the VICC JSON knowledgebase");
        options.addOption(ICLUSION_TRIAL_TSV, true, "Path to the iClusion input trial TSV");
        options.addOption(DOCM_TSV, true, "Path to the DoCM knowledgebase input TSV");
        options.addOption(HARTWIG_COHORT_TSV, true, "Path to the Hartwig Cohort input TSV");
        options.addOption(HARTWIG_CURATED_TSV, true, "Path to the Hartwig Curated input TSV");

        options.addOption(REF_GENOME_VERSION, true, "Ref version. Should be 'hgXX'");
        options.addOption(REF_GENOME_FASTA_FILE, true, "Path to the ref genome fasta file");
        options.addOption(DRIVER_GENE_TSV, true, "Path to driver gene TSV");

        options.addOption(OUTPUT_DIR, true, "Dir which will hold all SERVE output files");

        options.addOption(SKIP_HOTSPOT_RESOLVING, false, "If present, skips hotspot resolving");
        options.addOption(LOG_DEBUG, false, "If present, switches the logging to DEBUG mode");

        return options;
    }

    @NotNull
    String viccJson();

    @NotNull
    String iClusionTrialTsv();

    @NotNull
    String docmTsv();

    @NotNull
    String hartwigCohortTsv();

    @NotNull
    String hartwigCuratedTsv();

    @NotNull
    RefGenomeVersion refGenomeVersion();

    @NotNull
    String refGenomeFastaFile();

    @NotNull
    String driverGeneTsv();

    @NotNull
    String outputDir();

    boolean skipHotspotResolving();

    @NotNull
    static ServeConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        return ImmutableServeConfig.builder()
                .viccJson(nonOptionalFile(cmd, VICC_JSON))
                .iClusionTrialTsv(nonOptionalFile(cmd, ICLUSION_TRIAL_TSV))
                .docmTsv(nonOptionalFile(cmd, DOCM_TSV))
                .hartwigCohortTsv(nonOptionalFile(cmd, HARTWIG_COHORT_TSV))
                .hartwigCuratedTsv(nonOptionalFile(cmd, HARTWIG_CURATED_TSV))
                .refGenomeVersion(RefGenomeVersion.valueOf(nonOptionalValue(cmd, REF_GENOME_VERSION)))
                .refGenomeFastaFile(nonOptionalFile(cmd, REF_GENOME_FASTA_FILE))
                .driverGeneTsv(nonOptionalFile(cmd, DRIVER_GENE_TSV))
                .outputDir(nonOptionalDir(cmd, OUTPUT_DIR))
                .skipHotspotResolving(cmd.hasOption(SKIP_HOTSPOT_RESOLVING))
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
