package com.hartwig.hmftools.protect;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.protect.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.protect.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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

    private static final String OUTPUT_DATABASE_TSV = "output_database_tsv";
    private static final String OUTPUT_REPORT_TSV = "output_report_tsv";

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

        // Params output file
        final String outputDatabaseTsv = cmd.getOptionValue(OUTPUT_DATABASE_TSV);
        final String outputReportTsv = cmd.getOptionValue(OUTPUT_REPORT_TSV);

        if (!validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        LOGGER.info("Reading knowledgebase from {}", knowledgebaseDirectory);
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);

        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationCsv, tumorSampleId);
        List<? extends Variant> passSomaticVariants = readPassSomaticVariants(tumorSampleId, somaticVariantVcf);
        double ploidy = extractPloidy(purplePurityTsv);
        List<GeneCopyNumber> geneCopyNumbers = readGeneCopyNumbers(purpleGeneCnvTsv);
        List<ReportableGeneFusion> geneFusions = readGeneFusions(linxFusionTsv);

        LOGGER.info("Create actionability for sample {}", tumorSampleId);

        List<EvidenceItem> combinedEvidence = createEvidenceForAllFindings(actionabilityAnalyzer,
                patientPrimaryTumorLocation,
                passSomaticVariants,
                ploidy,
                geneCopyNumbers,
                geneFusions);

        LOGGER.info("Write actionability for patient report");
        writeActionabilityForPatientReport(outputReportTsv, combinedEvidence);

        LOGGER.info("Create actionability for database");
        writeActionabilityForDatabase(outputDatabaseTsv, combinedEvidence);

        LOGGER.info("Create conclusion for sample");
        writeConclusionOfSample();
    }

    private static void writeConclusionOfSample() {
        //TODO create conclusion
    }
    private static void writeActionabilityForPatientReport(@NotNull String outputReportTsv, @NotNull List<EvidenceItem> combinedEvidence)
            throws IOException {
        //TODO filter actionability

        BufferedWriter writerReport = new BufferedWriter(new FileWriter(outputReportTsv, false));
        writerReport.write(
                "event" + "\t" + "source" + "\t" + "reference" + "\t" + "drug" + "\t" + "drugsType" + "\t" + "level" + "\t" + "response"
                        + "\t" + "isOnLabel" + "\t" + "cancerType" + "\t" + "scope" + "\n");
        for (EvidenceItem item : combinedEvidence) {
            writerReport.write(
                    item.event() + "\t" + item.source() + "\t" + item.reference() + "\t" + item.drug() + "\t" + item.drugsType() + "\t"
                            + item.level() + "\t" + item.response() + "\t" + item.isOnLabel() + "\t" + item.cancerType() + "\t"
                            + item.scope() + "\n");
        }
        writerReport.close();
    }

    private static void writeActionabilityForDatabase(@NotNull String outputDatabaseTsv, @NotNull List<EvidenceItem> combinedEvidence)
            throws IOException {
        BufferedWriter writerDatabase = new BufferedWriter(new FileWriter(outputDatabaseTsv, false));
        writerDatabase.write(
                "event" + "\t" + "source" + "\t" + "reference" + "\t" + "drug" + "\t" + "drugsType" + "\t" + "level" + "\t" + "response"
                        + "\t" + "isOnLabel" + "\t" + "cancerType" + "\t" + "scope" + "\n");
        for (EvidenceItem item : combinedEvidence) {
            writerDatabase.write(
                    item.event() + "\t" + item.source() + "\t" + item.reference() + "\t" + item.drug() + "\t" + item.drugsType() + "\t"
                            + item.level() + "\t" + item.response() + "\t" + item.isOnLabel() + "\t" + item.cancerType() + "\t"
                            + item.scope() + "\n");
        }
        writerDatabase.close();
    }

    private static double extractPloidy(@NotNull String purplePurityTsv) throws IOException {
        LOGGER.info("Reading purple purity from {}", purplePurityTsv);
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info(" Sample ploidy: {}", ploidy);
        return ploidy;
    }

    @NotNull
    private static List<EvidenceItem> createEvidenceForAllFindings(@NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @NotNull String patientPrimaryTumorLocation, @NotNull List<? extends Variant> variants, double ploidy,
            @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull List<ReportableGeneFusion> geneFusions) {
        LOGGER.info("Extracting all evidence");

        List<EvidenceItem> evidenceForVariants =
                toList(actionabilityAnalyzer.evidenceForAllVariants(variants, patientPrimaryTumorLocation));
        LOGGER.info(" Found {} evidence items for {} somatic variants.", evidenceForVariants.size(), variants.size());

        List<EvidenceItem> evidenceForGeneCopyNumbers =
                toList(actionabilityAnalyzer.evidenceForCopyNumbers(geneCopyNumbers, patientPrimaryTumorLocation, ploidy));
        LOGGER.info(" Found {} evidence items for {} copy numbers.", evidenceForGeneCopyNumbers.size(), geneCopyNumbers.size());

        List<EvidenceItem> evidenceForGeneFusions =
                toList(actionabilityAnalyzer.evidenceForFusions(geneFusions, patientPrimaryTumorLocation));
        LOGGER.info(" Found {} evidence items for {} gene fusions.", evidenceForGeneFusions.size(), geneFusions.size());

        List<EvidenceItem> combinedEvidence = Lists.newArrayList();
        combinedEvidence.addAll(evidenceForVariants);
        combinedEvidence.addAll(evidenceForGeneCopyNumbers);
        combinedEvidence.addAll(evidenceForGeneFusions);
        return combinedEvidence;
    }

    @NotNull
    private static List<EvidenceItem> toList(@NotNull Map<?, List<EvidenceItem>> evidenceItemMap) {
        List<EvidenceItem> evidenceItemList = Lists.newArrayList();
        for (List<EvidenceItem> items : evidenceItemMap.values()) {
            evidenceItemList.addAll(items);
        }
        return evidenceItemList;
    }

    @NotNull
    private static String extractPatientTumorLocation(@NotNull String tumorLocationCsv, @NotNull String sampleId) throws IOException {
        LOGGER.info("Reading primary tumor location from {}", tumorLocationCsv);
        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info(" Loaded tumor locations for {} patients", patientTumorLocations.size());

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(patientTumorLocations, sampleId);

        String patientPrimaryTumorLocation = Strings.EMPTY;
        if (patientTumorLocation != null) {
            patientPrimaryTumorLocation = patientTumorLocation.primaryTumorLocation();
        }

        LOGGER.info(" Retrieved tumor location '{}' for sample {}", patientPrimaryTumorLocation, sampleId);

        return patientPrimaryTumorLocation;
    }

    @NotNull
    private static List<? extends Variant> readPassSomaticVariants(@NotNull String sampleId, @NotNull String somaticVariantVcf)
            throws IOException {
        LOGGER.info("Reading somatic variants from {}", somaticVariantVcf);
        List<? extends Variant> passSomaticVariants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sampleId, somaticVariantVcf);
        LOGGER.info(" Loaded {} PASS somatic variants", passSomaticVariants.size());
        return passSomaticVariants;
    }

    @NotNull
    private static List<GeneCopyNumber> readGeneCopyNumbers(@NotNull String purpleGeneCnvTsv) throws IOException {
        LOGGER.info("Reading gene copy numbers from {}", purpleGeneCnvTsv);
        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info(" Loaded {} gene copy numbers", geneCopyNumbers.size());
        return geneCopyNumbers;
    }

    @NotNull
    private static List<ReportableGeneFusion> readGeneFusions(@NotNull String linxFusionTsv) throws IOException {
        LOGGER.info("Reading gene fusions from {}", linxFusionTsv);
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info(" Loaded {} fusions", fusions.size());
        return fusions;
    }

    private static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, TUMOR_SAMPLE_ID) && dirExists(cmd, KNOWLEDGEBASE_DIRECTORY) && fileExists(cmd, TUMOR_LOCATION_CSV)
                && fileExists(cmd, SOMATIC_VARIANT_VCF) && fileExists(cmd, PURPLE_PURITY_TSV) && fileExists(cmd, PURPLE_GENE_CNV_TSV)
                && fileExists(cmd, LINX_FUSION_TSV);
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

        options.addOption(OUTPUT_DATABASE_TSV, true, "Path towards the output file for the database TSV.");
        options.addOption(OUTPUT_REPORT_TSV, true, "Path towards the output file for the report TSV.");

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
