package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.reportablegenomicalterations.AllReportableVariants;
import com.hartwig.hmftools.common.bachelor.BachelorFile;
import com.hartwig.hmftools.common.bachelor.FilterGermlineVariants;
import com.hartwig.hmftools.common.bachelor.GermlineReportingFile;
import com.hartwig.hmftools.common.bachelor.GermlineReportingModel;
import com.hartwig.hmftools.common.bachelor.GermlineVariant;
import com.hartwig.hmftools.common.reportablegenomicalterations.ReportableGainLoss;
import com.hartwig.hmftools.common.reportablegenomicalterations.ReportableGermlineVariant;
import com.hartwig.hmftools.common.reportablegenomicalterations.extractReportableGainsAndLosses;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.driverGene.DriverGeneView;
import com.hartwig.hmftools.common.driverGene.DriverGeneViewFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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

public class LoadEvidenceData {

    private static final Logger LOGGER = LogManager.getLogger(LoadEvidenceData.class);

    private static final String SAMPLE = "sample";
    private static final String TUMOR_SAMPLE_BARCODE = "tumor_sample_barcode";

    private static final String RUN_DIR = "run_dir";
    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String TUMOR_LOCATION_CSV = "tumorlocation_file";
    private static final String LIMS_DIRECTORY = "lims_dir";
    private static final String GERMLINE_GENES_CSV = "germline_genes_csv";

    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";
    private static final String BACHELOR_TSV = "bachelor_tsv";
    private static final String CHORD_PREDICTION_TXT = "chord_prediction_txt";
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        // Access for db
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);

        // General params needed for every sample
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        final String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);
        final String limsDirectory = cmd.getOptionValue(LIMS_DIRECTORY);
        final String germlineGenesCsv = cmd.getOptionValue(GERMLINE_GENES_CSV);

        // Params specific for specific sample
        final String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        final String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        final String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);
        final String bachelorTsv = cmd.getOptionValue(BACHELOR_TSV);
        final String chordPredictionTxt = cmd.getOptionValue(CHORD_PREDICTION_TXT);
        final String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final String tumorSampleBarcode = cmd.getOptionValue(TUMOR_SAMPLE_BARCODE);
        final String runDirectoryPath = cmd.getOptionValue(RUN_DIR);

        if (Utils.anyNull(userName,
                password,
                databaseUrl,
                runDirectoryPath,
                knowledgebaseDirectory,
                sampleId,
                tumorSampleBarcode,
                tumorLocationCsv,
                somaticVariantVcf,
                purpleGeneCnvTsv,
                linxFusionTsv,
                bachelorTsv,
                limsDirectory,
                chordPredictionTxt,
                germlineGenesCsv,
                purplePurityTsv)) {
            printUsageAndExit(options);
        }

        final File runDirectory = new File(runDirectoryPath);
        if (!runDirectory.isDirectory()) {
            LOGGER.warn("run_dir {} has to be an actual directory", runDirectory);
            printUsageAndExit(options);
        }

        LOGGER.info("Connecting with database");
        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Reading knowledgebase");
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);

        LOGGER.info("SampleId: " + sampleId);
        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationCsv, sampleId);

        List<ReportableGeneFusion> fusions = readingGeneFusions(linxFusionTsv);

        double ploidy = extractPloidy(purplePurityTsv);
        List<GeneCopyNumber> geneCopyNumbers = readingGeneCopyNumbers(purpleGeneCnvTsv);

        LOGGER.info("Extract reportable gains and losses ");
        List<ReportableGainLoss> reportableGainLosses = extractReportableGainsAndLosses.toReportableGainsAndLosses(geneCopyNumbers, ploidy);

        LOGGER.info("Extract all reportable variants of sample");
        List<SomaticVariant> passSomaticVariants = readingSomaticVariants(sampleId, somaticVariantVcf);
        List<GermlineVariant> passGermlineVariants = readingGermlineVariants(bachelorTsv);
        ChordAnalysis chordAnalysis = readingChord(chordPredictionTxt);

        LOGGER.info("Retrieve driver gene view");
        final DriverGeneView driverGeneView = DriverGeneViewFactory.create();

        LOGGER.info("Retrieve germline reporting model");
        final GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(germlineGenesCsv);

        List<ReportableVariant> extractAllReportableVariants = extarctReportableVariants(sampleId,
                limsDirectory,
                passSomaticVariants,
                passGermlineVariants,
                chordAnalysis,
                geneCopyNumbers,
                driverGeneView,
                germlineReportingModel,
                tumorSampleBarcode);

        List<EvidenceItem> combinedEvidence = createEvidenceOfAllFindings(actionabilityAnalyzer,
                patientPrimaryTumorLocation,
                extractAllReportableVariants,
                geneCopyNumbers,
                fusions,
                ploidy,
                reportableGainLosses);

        LOGGER.info("Writing evidence items into db");
        dbAccess.writeClinicalEvidence(sampleId, combinedEvidence);
        LOGGER.info("Finished");
    }

    private static double extractPloidy(@NotNull String purplePurityTsv) throws IOException {
        LOGGER.info("Reading purple purity file");
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info("Sample ploidy: " + ploidy);
        return ploidy;
    }

    @NotNull
    private static List<EvidenceItem> createEvidenceOfAllFindings(@NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @NotNull String patientPrimaryTumorLocation, @NotNull List<ReportableVariant> extractAllReportableVariants,
            @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull List<ReportableGeneFusion> fusions, double ploidy,
            @NotNull List<ReportableGainLoss> reportableGainLosses) {
        LOGGER.info("Exctracting all evidence");
        Map<ReportableVariant, List<EvidenceItem>> evidencePerVariant =
                actionabilityAnalyzer.evidenceForAllVariants(extractAllReportableVariants, patientPrimaryTumorLocation);

        List<EvidenceItem> allEvidenceForVariants = extractAllEvidenceItems(evidencePerVariant);
        LOGGER.info("Found {} evidence items for {} somatic variants.", allEvidenceForVariants.size(), extractAllReportableVariants.size());

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(geneCopyNumbers, patientPrimaryTumorLocation, ploidy);
        List<EvidenceItem> allEvidenceForCopyNumbers = extractAllEvidenceItems(evidencePerGeneCopyNumber);

        //TODO fix is semi duplicate code as is CopyNumberAnalysis!
        // Check that all copy numbers with evidence are reported (since they are in the driver catalog).
        Set<String> reportableGenes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            reportableGenes.add(gainLoss.gene());
        }

        for (Map.Entry<GeneCopyNumber, List<EvidenceItem>> entry : evidencePerGeneCopyNumber.entrySet()) {
            GeneCopyNumber geneCopyNumber = entry.getKey();
            if (!Collections.disjoint(entry.getValue(), allEvidenceForCopyNumbers) && !reportableGenes.contains(geneCopyNumber.gene())) {
                LOGGER.warn("Copy number with evidence not reported: {}!", geneCopyNumber.gene());
            }
        }

        LOGGER.info("Found {} evidence items for {} copy numbers.",
                allEvidenceForCopyNumbers.size(),
                evidencePerGeneCopyNumber.keySet().size());

        Map<ReportableGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(fusions, patientPrimaryTumorLocation);
        List<EvidenceItem> allEvidenceForGeneFusions = extractAllEvidenceItems(evidencePerFusion);
        LOGGER.info("Found {} evidence items for {} gene fusions.", allEvidenceForGeneFusions.size(), fusions.size());

        List<EvidenceItem> combinedEvidence = Lists.newArrayList();
        combinedEvidence.addAll(allEvidenceForVariants);
        combinedEvidence.addAll(allEvidenceForCopyNumbers);
        combinedEvidence.addAll(allEvidenceForGeneFusions);
        return combinedEvidence;
    }

    @NotNull
    private static ChordAnalysis readingChord(@NotNull String chordPredictionTxt) throws IOException {
        LOGGER.info("Reading chord from file");
        ChordAnalysis chordAnalysis = ChordFileReader.read(chordPredictionTxt);
        LOGGER.info("Loaded CHORD analysis from {}", chordPredictionTxt);
        return chordAnalysis;
    }

    @NotNull
    private static List<GermlineVariant> readingGermlineVariants(@NotNull String bachelorTsv) throws IOException {
        LOGGER.info("Reading germline variants from file");
        List<GermlineVariant> passGermlineVariants =
                BachelorFile.loadBachelorTsv(bachelorTsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", passGermlineVariants.size(), bachelorTsv);
        return passGermlineVariants;
    }

    @NotNull
    private static List<SomaticVariant> readingSomaticVariants(@NotNull String sampleId, @NotNull String somaticVariantVcf)
            throws IOException {
        LOGGER.info("Reading somatic variants from file");
        List<SomaticVariant> passSomaticVariants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sampleId, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", passSomaticVariants.size(), somaticVariantVcf);
        return passSomaticVariants;
    }

    @NotNull
    private static List<ReportableVariant> extarctReportableVariants(@NotNull String sampleId, @NotNull String limsDirectory,
            @NotNull List<SomaticVariant> passSomaticVariants, @NotNull List<GermlineVariant> passGermlineVariants,
            @NotNull ChordAnalysis chordAnalysis, List<GeneCopyNumber> geneCopyNumbers, @NotNull DriverGeneView driverGeneView,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull String tumorSampleBarcode) throws IOException {

        Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);

        LOGGER.info("Filtering germline variants");
        List<ReportableGermlineVariant> reportableGermlineVariants = filterGermlineVariants(sampleId,
                limsDirectory,
                chordAnalysis,
                passSomaticVariants,
                passGermlineVariants,
                geneCopyNumbers,
                driverGeneView,
                germlineReportingModel,
                lims);

        List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(passSomaticVariants, geneCopyNumbers));
        driverCatalog.addAll(TsgDrivers.drivers(passSomaticVariants, geneCopyNumbers));

        LOGGER.info("Merging all reportable somatic and germline variants");
        return AllReportableVariants.mergeSomaticAndGermlineVariants(passSomaticVariants,
                driverCatalog,
                driverGeneView,
                reportableGermlineVariants,
                germlineReportingModel,
                lims.germlineReportingChoice(tumorSampleBarcode));
    }

    @NotNull
    private static List<ReportableGermlineVariant> filterGermlineVariants(@NotNull String sampleId, @NotNull String limsDirectory,
            @NotNull ChordAnalysis chordAnalysis, @NotNull List<SomaticVariant> passSomaticVariants,
            @NotNull List<GermlineVariant> passGermlineVariants, List<GeneCopyNumber> geneCopyNumbers,
            @NotNull DriverGeneView driverGeneView, @NotNull GermlineReportingModel germlineReportingModel, @NotNull Lims lims)
            throws IOException {

        LimsGermlineReportingChoice germlineReportingChoice = extarctGermlineChoiceOfPatient(sampleId, limsDirectory, lims);
        if (germlineReportingChoice == LimsGermlineReportingChoice.UNKNOWN) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
            return Lists.newArrayList();
        } else {
            LOGGER.info(" Patient has given the following germline consent: {}", germlineReportingChoice);
            return FilterGermlineVariants.filterGermlineVariantsForReporting(passGermlineVariants,
                    driverGeneView,
                    germlineReportingModel,
                    geneCopyNumbers,
                    passSomaticVariants,
                    chordAnalysis);
        }
    }

    @NotNull
    private static List<GeneCopyNumber> readingGeneCopyNumbers(@NotNull String purpleGeneCnvTsv) throws IOException {
        LOGGER.info("Reading gene copy numbers from file");
        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info("Loaded {} gene copy numbers from {}", geneCopyNumbers.size(), purpleGeneCnvTsv);
        return geneCopyNumbers;
    }

    @NotNull
    private static List<ReportableGeneFusion> readingGeneFusions(@NotNull String linxFusionTsv) throws IOException {
        LOGGER.info("Reading gene fusions from file");
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.read(linxFusionTsv);
        LOGGER.info("Loaded {} fusions from {}", fusions.size(), linxFusionTsv);
        return fusions;
    }

    @NotNull
    private static String extractPatientTumorLocation(@NotNull String tumorLocationCsv, @NotNull String sampleId) throws IOException {
        LOGGER.info("Reading primary tumor location from file.");
        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(tumorLocationCsv);
        LOGGER.info("Loaded tumor locations for {} patients from {}", patientTumorLocations.size(), tumorLocationCsv);

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(patientTumorLocations, sampleId);

        String patientPrimaryTumorLocation = Strings.EMPTY;
        if (patientTumorLocation != null) {
            patientPrimaryTumorLocation = patientTumorLocation.primaryTumorLocation();
        }

        LOGGER.info("Primary tumor location of patient: " + patientPrimaryTumorLocation);
        return patientPrimaryTumorLocation;
    }

    @NotNull
    private static LimsGermlineReportingChoice extarctGermlineChoiceOfPatient(@NotNull String sampleId, @NotNull String limsDirectory,
            @NotNull Lims lims) throws IOException {
        LOGGER.info("Loaded LIMS data for {} samples from {}", lims.sampleBarcodeCount(), limsDirectory);
        LimsGermlineReportingChoice germlineChoice = lims.germlineReportingChoice(sampleId);
        if (germlineChoice == LimsGermlineReportingChoice.UNKNOWN) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
        } else {
            LOGGER.info(" Patient has given the following germline consent: {}", germlineChoice);
        }
        return germlineChoice;
    }

    @NotNull
    private static List<EvidenceItem> extractAllEvidenceItems(@NotNull Map<?, List<EvidenceItem>> evidenceItemMap) {
        List<EvidenceItem> evidenceItemList = Lists.newArrayList();
        for (List<EvidenceItem> items : evidenceItemMap.values()) {
            evidenceItemList.addAll(items);
        }
        return evidenceItemList;
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("patient-db - load evidence data", options);
        System.exit(1);
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the folder containing knowledgebase files.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");
        options.addOption(LIMS_DIRECTORY, true, "Path towards the directory holding the LIMS data");
        options.addOption(GERMLINE_GENES_CSV, true, "Path towards a CSV containing germline genes which we want to report.");

        options.addOption(SAMPLE, true, "Tumor sample of run");
        options.addOption(TUMOR_SAMPLE_BARCODE, true, "The sample barcode for which a patient report will be generated.");
        options.addOption(RUN_DIR, true, "Path towards the folder containing sample run.");
        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(BACHELOR_TSV, true, "Path towards the germline TSV (optional).");
        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT .");
        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
