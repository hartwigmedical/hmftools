package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.bachelor.BachelorFile;
import com.hartwig.hmftools.common.bachelor.GermlineVariant;
import com.hartwig.hmftools.common.bachelor.ReportableGermlineVariant;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.lims.LimsFactory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Drivercatalog;

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
        final String sampleId = cmd.getOptionValue(SAMPLE);
        final String runDirectoryPath = cmd.getOptionValue(RUN_DIR);

        if (Utils.anyNull(userName,
                password,
                databaseUrl,
                runDirectoryPath,
                knowledgebaseDirectory,
                sampleId,
                tumorLocationCsv,
                somaticVariantVcf,
                purpleGeneCnvTsv,
                linxFusionTsv,
                bachelorTsv,
                limsDirectory,
                chordPredictionTxt,
                germlineGenesCsv)) {
            printUsageAndExit(options);
        }

        final File runDirectory = new File(runDirectoryPath);
        if (!runDirectory.isDirectory()) {
            LOGGER.warn("run_dir {} has to be an actual directory", runDirectory);
            printUsageAndExit(options);
        }

        DatabaseAccess dbAccess = databaseAccess(cmd);

        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);
        //  final GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(germlineGenesCsv);
        //        final DriverGeneView driverGeneView = DriverGeneViewFactory.create();
        //
        List<Drivercatalog> driverCatalog = Lists.newArrayList();
        //        driverCatalog.addAll(OncoDrivers.drivers(variants, exomeGeneCopyNumbers));
        //        driverCatalog.addAll(TsgDrivers.drivers(variants, exomeGeneCopyNumbers));

        LOGGER.info("Sample: " + sampleId);

        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationCsv, sampleId);

        LOGGER.info("Reading gene copy numbers and sample ploidy from DB");

        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info("Loaded {} gene copy numbers from {}", geneCopyNumbers.size(), purpleGeneCnvTsv);
        List<ReportableGeneFusion> fusions = readingGeneFusions(linxFusionTsv);
        List<SomaticVariant> passSomaticVariants =
                extarctReportableVariants(sampleId, limsDirectory, somaticVariantVcf, bachelorTsv, chordPredictionTxt, geneCopyNumbers);

        Map<SomaticVariant, List<EvidenceItem>> evidencePerVariant =
                actionabilityAnalyzer.evidenceForSomaticVariants(passSomaticVariants, patientPrimaryTumorLocation);

        List<EvidenceItem> allEvidenceForSomaticVariants = extractAllEvidenceItems(evidencePerVariant);
        LOGGER.info("Found {} evidence items for {} somatic variants.", allEvidenceForSomaticVariants.size(), passSomaticVariants.size());

        PurityContext purityContext = dbAccess.readPurityContext(sampleId);
        assert purityContext != null;

        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info("Sample ploidy: " + ploidy);

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(geneCopyNumbers, patientPrimaryTumorLocation, ploidy);
        List<EvidenceItem> allEvidenceForCopyNumbers = extractAllEvidenceItems(evidencePerGeneCopyNumber);
        LOGGER.info("Found {} evidence items for {} copy numbers.",
                allEvidenceForCopyNumbers.size(),
                evidencePerGeneCopyNumber.keySet().size());

        Map<ReportableGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(fusions, patientPrimaryTumorLocation);
        List<EvidenceItem> allEvidenceForGeneFusions = extractAllEvidenceItems(evidencePerFusion);
        LOGGER.info("Found {} evidence items for {} gene fusions.", allEvidenceForGeneFusions.size(), fusions.size());

        List<EvidenceItem> combinedEvidence = Lists.newArrayList();
        combinedEvidence.addAll(allEvidenceForSomaticVariants);
        combinedEvidence.addAll(allEvidenceForCopyNumbers);
        combinedEvidence.addAll(allEvidenceForGeneFusions);

        dbAccess.writeClinicalEvidence(sampleId, combinedEvidence);
    }

    @NotNull
    private static List<SomaticVariant> extarctReportableVariants(@NotNull String sampleId, @NotNull String limsDirectory,
            @NotNull String somaticVariantVcf, @NotNull String bachelorTsv, @NotNull String chordPredictionTxt,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) throws IOException {

        LOGGER.info("Reading somatic variants from DB");
        List<SomaticVariant> passSomaticVariants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sampleId, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", passSomaticVariants.size(), somaticVariantVcf);
        List<ReportableGermlineVariant> reportableGermlineVariants =
                filterGermlineVariants(bachelorTsv, sampleId, limsDirectory, chordPredictionTxt, passSomaticVariants, geneCopyNumbers);
        // merge somatic and germline variants
        return passSomaticVariants;
    }

    @NotNull
    private static List<ReportableGermlineVariant> filterGermlineVariants(@NotNull String bachelorTsv, @NotNull String sampleId,
            @NotNull String limsDirectory, @NotNull String chordPredictionTxt, @NotNull List<SomaticVariant> passSomaticVariants,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) throws IOException {
        ChordAnalysis chordAnalysis = ChordFileReader.read(chordPredictionTxt);
        LOGGER.info("Loaded CHORD analysis from {}", chordPredictionTxt);

        LOGGER.info("Reading germline variants from DB");
        List<GermlineVariant> germlineVariants =
                BachelorFile.loadBachelorTsv(bachelorTsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info("Loaded {} PASS germline variants from {}", germlineVariants.size(), bachelorTsv);

        LimsGermlineReportingChoice germlineReportingChoice = extarctGermlineChoiceOfPatient(sampleId, limsDirectory);
        if (germlineReportingChoice == LimsGermlineReportingChoice.UNKNOWN) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
            return Lists.newArrayList();
        } else {
            LOGGER.info(" Patient has given the following germline consent: {}", germlineReportingChoice);
            return Lists.newArrayList();
            //            return FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariants,
            //                    reportData.driverGeneView(),
            //                    reportData.germlineReportingModel(),
            //                    geneCopyNumbers,
            //                    passSomaticVariants,
            //                    chordAnalysis);
        }

    }

    @NotNull
    private static List<ReportableGeneFusion> readingGeneFusions(@NotNull String linxFusionTsv) throws IOException {
        LOGGER.info("Reading gene fusions from DB");
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

        LOGGER.info("Primary tumor location of patient: " + patientTumorLocation);
        return patientPrimaryTumorLocation;
    }

    @NotNull
    private static LimsGermlineReportingChoice extarctGermlineChoiceOfPatient(@NotNull String sampleId, @NotNull String limsDirectory)
            throws IOException {
        Lims lims = LimsFactory.fromLimsDirectory(limsDirectory);
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
        options.addOption(RUN_DIR, true, "Path towards the folder containing sample run.");
        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");
        options.addOption(BACHELOR_TSV, true, "Path towards the germline TSV (optional).");
        options.addOption(CHORD_PREDICTION_TXT, true, "Path towards the CHORD prediction TXT .");

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
