package com.hartwig.hmftools.protect;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalLoad;
import com.hartwig.hmftools.protect.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.protect.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.protect.common.GenomicData;
import com.hartwig.hmftools.protect.conclusion.ConclusionFactory;
import com.hartwig.hmftools.protect.conclusion.TemplateConclusion;
import com.hartwig.hmftools.protect.conclusion.TemplateConclusionFile;
import com.hartwig.hmftools.protect.report.chord.ChordFileReader;

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

public class ProtectActionability {
    private static final Logger LOGGER = LogManager.getLogger(ProtectActionability.class);

    private static final String TUMOR_SAMPLE_ID = "tumor_sample_id";

    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String TEMPLATE_CONCLUSION_TSV = "template_conclusion";

    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String CHORD_TXT = "chord_txt";

    private static final String CONCLUSION_TSV = "conclusion_tsv";
    private static final String OUTPUT_DATABASE_TSV = "output_database_tsv";
    private static final String OUTPUT_REPORT_TSV = "output_report_tsv";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        String tumorSampleId = cmd.getOptionValue(TUMOR_SAMPLE_ID);

        // General params needed for every sample
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        final String templateConclusionTsv = cmd.getOptionValue(TEMPLATE_CONCLUSION_TSV);

        // Params specific for specific sample
        final String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        final String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);
        final String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        final String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);
        final String chordTxt = cmd.getOptionValue(CHORD_TXT);

        // Params output file
        final String outputDatabaseTsv = cmd.getOptionValue(OUTPUT_DATABASE_TSV);
        final String outputReportTsv = cmd.getOptionValue(OUTPUT_REPORT_TSV);
        final String OutputConclusionTsv = cmd.getOptionValue(CONCLUSION_TSV);

        if (!validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        LOGGER.info("Reading knowledgebase from {}", knowledgebaseDirectory);
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);

        LOGGER.info("Reading template Conclusion from {}", templateConclusionTsv);
        List<TemplateConclusion> templateConclusionList = TemplateConclusionFile.readTemplateConclusion(templateConclusionTsv);

        // Extract genomic alterations
        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationCsv, tumorSampleId);
        List<? extends Variant> passSomaticVariants = GenomicData.readPassSomaticVariants(tumorSampleId, somaticVariantVcf);
        double ploidy = GenomicData.extractPloidy(purplePurityTsv);
        List<GeneCopyNumber> geneCopyNumbers = GenomicData.readGeneCopyNumbers(purpleGeneCnvTsv);
        List<ReportableGeneFusion> geneFusions = GenomicData.readGeneFusions(linxFusionTsv);

        // Extract tumor characteristics
        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(tumorSampleId, somaticVariantVcf);
        double tumorMTL = TumorMutationalLoad.determineTumorMutationalLoad(variants);
        double tumorMSI = MicrosatelliteIndels.determineMicrosatelliteIndelsPerMb(variants);
        double chordScore = ChordFileReader.read(chordTxt).hrdValue();
        double tumorMTB = TumorMutationalLoad.determineTumorMutationalBurden(variants);

        LOGGER.info("Create actionability for sample {}", tumorSampleId);

        List<EvidenceItem> combinedEvidence = createEvidenceForAllFindings(actionabilityAnalyzer,
                patientPrimaryTumorLocation,
                passSomaticVariants,
                ploidy,
                geneCopyNumbers,
                geneFusions);

        LOGGER.info("Create actionability for patient report");
        writeActionabilityForPatientReport(outputReportTsv, combinedEvidence);

        LOGGER.info("Create actionability for database");
        writeActionabilityForDatabase(outputDatabaseTsv, combinedEvidence);

        LOGGER.info("Create conclusion for sample");
        String conclusion = ConclusionFactory.createConclusion(patientPrimaryTumorLocation,
                tumorMTL,
                tumorMTB,
                tumorMSI,
                chordScore,
                geneFusions,
                geneCopyNumbers, passSomaticVariants, templateConclusionList);

        writeConclusionOfSample(OutputConclusionTsv, conclusion);
    }

    private static void writeConclusionOfSample(@NotNull String OutputConclusionTsv, @NotNull String conclusion) throws IOException {
        //TODO create conclusion
        BufferedWriter writerReport = new BufferedWriter(new FileWriter(OutputConclusionTsv, false));

        writerReport.write(conclusion);

        writerReport.close();
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

    private static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, TUMOR_SAMPLE_ID) && dirExists(cmd, KNOWLEDGEBASE_DIRECTORY) && fileExists(cmd, TUMOR_LOCATION_CSV)
                && fileExists(cmd, SOMATIC_VARIANT_VCF) && fileExists(cmd, PURPLE_PURITY_TSV) && fileExists(cmd, PURPLE_GENE_CNV_TSV)
                && fileExists(cmd, LINX_FUSION_TSV) && fileExists(cmd, CHORD_TXT) && fileExists(cmd, TEMPLATE_CONCLUSION_TSV);
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
        options.addOption(TEMPLATE_CONCLUSION_TSV, true, "Path towards the template for conclusion TSV.");

        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(CHORD_TXT, true, "Path towards the chord txt file.");

        options.addOption(CONCLUSION_TSV, true, "Path towards the conclusion TSV.");
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
