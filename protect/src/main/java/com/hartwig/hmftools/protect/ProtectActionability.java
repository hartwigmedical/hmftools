package com.hartwig.hmftools.protect;

import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.cli.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.purple.copynumber.ExtractReportableGainsAndLosses;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.protect.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.protect.actionability.EvidenceItem;
import com.hartwig.hmftools.protect.common.BachelorFile;
import com.hartwig.hmftools.protect.common.CopyNumberAnalysis;
import com.hartwig.hmftools.protect.common.CopyNumberAnalyzer;
import com.hartwig.hmftools.protect.common.GenomicData;
import com.hartwig.hmftools.protect.common.GermlineReportingFile;
import com.hartwig.hmftools.protect.common.GermlineReportingModel;
import com.hartwig.hmftools.protect.common.GermlineVariant;
import com.hartwig.hmftools.protect.common.HomozygousDisruptionAnalyzer;
import com.hartwig.hmftools.protect.common.ReportableGermlineVariant;
import com.hartwig.hmftools.protect.common.ReportableHomozygousDisruption;
import com.hartwig.hmftools.protect.common.ReportableVariantAnalysis;
import com.hartwig.hmftools.protect.common.SomaticVariantAnalysis;
import com.hartwig.hmftools.protect.common.SomaticVariantAnalyzer;
import com.hartwig.hmftools.protect.conclusion.ConclusionFactory;
import com.hartwig.hmftools.protect.conclusion.TemplateConclusion;
import com.hartwig.hmftools.protect.conclusion.TemplateConclusionFile;

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
import org.jetbrains.annotations.Nullable;

public class ProtectActionability {

    private static final Logger LOGGER = LogManager.getLogger(ProtectActionability.class);

    private static final String TUMOR_SAMPLE_ID = "tumor_sample_id";
    private static final String TUMOR_BARCODE_ID = "tumor_barcode_id";

    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String KNOWLEDGEBASE_DIRECTORY_V2 = "knowledgebase_dir_v2";

    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
    private static final String TEMPLATE_CONCLUSION_TSV = "template_conclusion";
    private static final String GERMLINE_GENES_CSV = "germline_genes_csv";
    private static final String LIMS_DIRECTORY = "lims_dir";

    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String GERMLINE_VARIANT_VCF = "germline_variant_vcf";

    private static final String PURPLE_DRIVER_CATALOG_TSV = "purple_driver_catalog_tsv";
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String PURPLE_QC_TSV = "purple_qc_tsv";
    private static final String LINX_DRIVERS_TSV = "linx_drivers_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String CHORD_TXT = "chord_txt";

    private static final String CONCLUSION_TSV = "conclusion_tsv";
    private static final String OUTPUT_DATABASE_TSV = "output_database_tsv";
    private static final String OUTPUT_REPORT_TSV = "output_report_tsv";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        String tumorSampleId = cmd.getOptionValue(TUMOR_SAMPLE_ID);
        String tumorBarcodeId = cmd.getOptionValue(TUMOR_BARCODE_ID);

        // General params needed for every sample
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        final String knowledgebaseDirectory_v2 = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY_V2);

        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        final String templateConclusionTsv = cmd.getOptionValue(TEMPLATE_CONCLUSION_TSV);
        final String germlineGenesCsv = cmd.getOptionValue(GERMLINE_GENES_CSV);
        final String limsDir = cmd.getOptionValue(LIMS_DIRECTORY);

        // Params specific for specific sample
        final String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        final String germlineVariantVcf = cmd.getOptionValue(GERMLINE_VARIANT_VCF);
        final String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);
        final String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        final String purpleQCTsv = cmd.getOptionValue(PURPLE_QC_TSV);
        final String linxDriversTsv = cmd.getOptionValue(LINX_DRIVERS_TSV);
        final String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);
        final String chordTxt = cmd.getOptionValue(CHORD_TXT);

        final String driverGenePanelFile = cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION);

        // Params output file
        final String outputDatabaseTsv = cmd.getOptionValue(OUTPUT_DATABASE_TSV);
        final String outputReportTsv = cmd.getOptionValue(OUTPUT_REPORT_TSV);
        final String OutputConclusionTsv = cmd.getOptionValue(CONCLUSION_TSV);

        if (!validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        LOGGER.info("Reading knowledgebase generator v2 from {}", knowledgebaseDirectory_v2);
        com.hartwig.hmftools.serve.actionability.ActionabilityAnalyzer actionabilityAnalyzer_v2 =
                com.hartwig.hmftools.serve.actionability.ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory_v2);

        LOGGER.info("Reading knowledgebase from {}", knowledgebaseDirectory);
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);

        LOGGER.info("Reading template Conclusion from {}", templateConclusionTsv);
        List<TemplateConclusion> templateConclusionList = TemplateConclusionFile.readTemplateConclusion(templateConclusionTsv);
        Map<String, TemplateConclusion> mapFindingToConclusion = Maps.newHashMap();

        for (TemplateConclusion templateConclusion : templateConclusionList) {
            mapFindingToConclusion.put(templateConclusion.abberrationGeneSummary(), templateConclusion);
        }

        LOGGER.info("Loading sample data from LIMS in {}", cmd.getOptionValue(LIMS_DIRECTORY));
        Lims lims = LimsFactory.fromLimsDirectory(cmd.getOptionValue(LIMS_DIRECTORY));

        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationCsv, tumorSampleId);
        String patientCancerSubtype = extractCancerSubtype(tumorLocationCsv, tumorSampleId);

        LOGGER.info("Extracting genomic alteration from sample {}", tumorSampleId);
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        double purity = purityContext.bestFit().purity();
        double ploidy = GenomicData.extractPloidy(purplePurityTsv);

        List<DriverGene> driverGenes = DriverGeneFile.read(driverGenePanelFile);
        List<DriverCatalog> driverCatalog = readDriverCatalog(cmd.getOptionValue(PURPLE_DRIVER_CATALOG_TSV));

        // Gene Fusion reportable
        List<ReportableGeneFusion> geneFusions = GenomicData.readGeneFusions(linxFusionTsv);

        // Copy Number all + reportable
        List<GeneCopyNumber> geneCopyNumbers = GenomicData.readGeneCopyNumbers(purpleGeneCnvTsv);
        List<ReportableGainLoss> reportableGainsAndLosses = ExtractReportableGainsAndLosses.toReportableGainsAndLosses(driverCatalog,
                geneCopyNumbers);

        // Germline variants
        List<GermlineVariant> germlineVariant =
                BachelorFile.loadBachelorTsv(germlineVariantVcf).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", germlineVariant.size(), germlineVariantVcf);

        // only reportable variants
        CopyNumberAnalysis copyNumberAnalysis =
                CopyNumberAnalyzer.analyzeCopyNumbers(purplePurityTsv, purpleQCTsv, purpleGeneCnvTsv, driverCatalog);
        //

        SomaticVariantAnalysis somaticVariantAnalysis = SomaticVariantAnalyzer.analyzeSomaticVariants(tumorSampleId,
                somaticVariantVcf,
                driverCatalog);
        ChordAnalysis chordAnalysis = ChordFileReader.read(chordTxt);
        final GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(germlineGenesCsv);

        List<ReportableGermlineVariant> germlineVariantsToReport = GenomicData.analyzeGermlineVariants(tumorBarcodeId,
                germlineVariantVcf,
                copyNumberAnalysis,
                somaticVariantAnalysis,
                chordAnalysis,
                driverGenes,
                lims,
                germlineReportingModel);

        ReportableVariantAnalysis reportableVariantsAnalysis =
                GenomicData.somaticAndGermlineVariantsTogether(somaticVariantAnalysis.variantsToReport(),
                        somaticVariantAnalysis.driverCatalog(),
                        driverGenes,
                        germlineVariantsToReport,
                        germlineReportingModel,
                        lims.germlineReportingChoice(tumorBarcodeId));

        List<ReportableHomozygousDisruption> reportableHomozygousDisruptions =
                HomozygousDisruptionAnalyzer.extractFromLinxDriversTsv(linxDriversTsv);

        LOGGER.info("Extract tumor characteristics from sample {}", tumorSampleId);
        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(tumorSampleId, somaticVariantVcf);
        int tumorMTL = (int) Math.round(copyNumberAnalysis.tumorMutationalLoad());
        double tumorMSI = copyNumberAnalysis.microsatelliteIndelsPerMb();
        double tumorMTB = copyNumberAnalysis.tumorMutationalBurdenPerMb();

        LOGGER.info("Create actionability for sample {}", tumorSampleId);

        List<EvidenceItem> combinedEvidence = Lists.newArrayList();
        //        List<EvidenceItem> combinedEvidence = createEvidenceForAllFindings(actionabilityAnalyzer,
        //                patientPrimaryTumorLocation,
        //                passSomaticVariants,
        //                ploidy,
        //                geneCopyNumbers,
        //                geneFusions);

        LOGGER.info("Create actionability for patient report");
        writeActionabilityForPatientReport(outputReportTsv, combinedEvidence);

        LOGGER.info("Create actionability for database");
        writeActionabilityForDatabase(outputDatabaseTsv, combinedEvidence);

        LOGGER.info("Create conclusion for sample");
        StringBuilder conclusion = ConclusionFactory.createConclusion(patientPrimaryTumorLocation,
                tumorMTL,
                tumorMTB,
                tumorMSI,
                chordAnalysis.hrdValue(),
                geneFusions,
                reportableGainsAndLosses,
                reportableVariantsAnalysis,
                purity,
                patientCancerSubtype,
                reportableHomozygousDisruptions,
                mapFindingToConclusion);

        LOGGER.info("Create hotspot information");
        //TODO create hotspot information

        LOGGER.info("Create knwon amplications and deletions");
        //TODO create knwon amplifications and deletions

        writeConclusionOfSample(OutputConclusionTsv, conclusion);
    }

    @NotNull
    public static List<DriverCatalog> readDriverCatalog(@NotNull String purpleDriverCatalogTsv) throws IOException {
        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(purpleDriverCatalogTsv);
        LOGGER.info("Loaded {} driver catalog records", driverCatalog.size());
        return driverCatalog;
    }

    private static void writeConclusionOfSample(@NotNull String OutputConclusionTsv, @NotNull StringBuilder conclusion) throws IOException {
        //TODO create conclusion
        BufferedWriter writerReport = new BufferedWriter(new FileWriter(OutputConclusionTsv, false));

        writerReport.write(conclusion.toString());

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
                            + item.level() + "\t" + item.response() + "\t" + item.isOnLabel() + "\t" + item.cancerType() + "\t" + "\n");
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
                            + item.level() + "\t" + item.response() + "\t" + item.isOnLabel() + "\t" + item.cancerType() + "\t" + "\n");
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

    @Nullable
    private static PatientTumorLocation extractTumorLocation(@NotNull String tumorLocationCsv, @NotNull String sampleId)
            throws IOException {
        LOGGER.info("Reading primary tumor location from {}", tumorLocationCsv);
        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info(" Loaded tumor locations for {} patients", patientTumorLocations.size());

        return PatientTumorLocationFunctions.findPatientTumorLocationForSample(patientTumorLocations, sampleId);
    }

    @NotNull
    private static String extractCancerSubtype(@NotNull String tumorLocationCsv, @NotNull String sampleId) throws IOException {
        PatientTumorLocation tumorLocation = extractTumorLocation(tumorLocationCsv, sampleId);

        String cancerSubtype = Strings.EMPTY;
        if (tumorLocation != null) {
            cancerSubtype = tumorLocation.cancerSubtype();
        }

        LOGGER.info(" Retrieved cancer subtupe '{}' for sample {}", cancerSubtype, sampleId);

        return cancerSubtype;
    }

    @NotNull
    private static String extractPatientTumorLocation(@NotNull String tumorLocationCsv, @NotNull String sampleId) throws IOException {
        PatientTumorLocation tumorLocation = extractTumorLocation(tumorLocationCsv, sampleId);

        String patientPrimaryTumorLocation = Strings.EMPTY;
        if (tumorLocation != null) {
            patientPrimaryTumorLocation = tumorLocation.primaryTumorLocation();
        }

        LOGGER.info(" Retrieved tumor location '{}' for sample {}", patientPrimaryTumorLocation, sampleId);

        return patientPrimaryTumorLocation;
    }

    private static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, TUMOR_SAMPLE_ID) && valueExists(cmd, TUMOR_BARCODE_ID) && dirExists(cmd, KNOWLEDGEBASE_DIRECTORY)
                && dirExists(cmd, KNOWLEDGEBASE_DIRECTORY_V2) && dirExists(cmd, LIMS_DIRECTORY) && fileExists(cmd, TUMOR_LOCATION_CSV)
                && fileExists(cmd, SOMATIC_VARIANT_VCF) && fileExists(cmd, PURPLE_PURITY_TSV) && fileExists(cmd, PURPLE_GENE_CNV_TSV)
                && fileExists(cmd, LINX_FUSION_TSV) && fileExists(cmd, CHORD_TXT) && fileExists(cmd, TEMPLATE_CONCLUSION_TSV) && fileExists(
                cmd,
                GERMLINE_VARIANT_VCF) && fileExists(cmd, PURPLE_QC_TSV) && fileExists(cmd, GERMLINE_GENES_CSV) && fileExists(cmd,
                LINX_DRIVERS_TSV);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            LOGGER.warn("'{}' has to be provided", param);
            return false;
        }
        return true;
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn("'{}' has to be an existing file: {}", param, value);
            return false;
        }

        return true;
    }

    private static boolean dirExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value) || !pathIsDirectory(value)) {
            LOGGER.warn("'{}' has to be an existing directory: {}", param, value);
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

        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);

        options.addOption(TUMOR_SAMPLE_ID, true, "The sample ID for which a patient report will be generated.");
        options.addOption(TUMOR_BARCODE_ID, true, "The barcode ID for which a patient report will be generated.");

        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the folder containing knowledgebase files.");
        options.addOption(KNOWLEDGEBASE_DIRECTORY_V2, true, "Path towards the folder containing knowledgebase files from version 2.");

        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");
        options.addOption(TEMPLATE_CONCLUSION_TSV, true, "Path towards the template for conclusion TSV.");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the LIMS directory.");

        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(GERMLINE_VARIANT_VCF, true, "Path towards the germline variant VCF.");

        options.addOption(PURPLE_DRIVER_CATALOG_TSV, true, "Path towards the purple driver catalog TSV.");
        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_TSV, true, "Path towards the purple qc TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(LINX_DRIVERS_TSV, true, "Path towards the LINX driver catalog TSV.");
        options.addOption(CHORD_TXT, true, "Path towards the chord txt file.");

        options.addOption(CONCLUSION_TSV, true, "Path towards the conclusion TSV.");
        options.addOption(OUTPUT_DATABASE_TSV, true, "Path towards the output file for the database TSV.");
        options.addOption(OUTPUT_REPORT_TSV, true, "Path towards the output file for the report TSV.");

        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        return new Options();
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Protect", options);
        System.exit(1);
    }
}
