package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusion;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusionFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
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

    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";
    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";

    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_GENE_CNV_TSV = "purple_gene_cnv_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(args, options);

        String sampleId = cmd.getOptionValue(SAMPLE);

        // General params needed for every sample
        String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        String tumorLocationCsv = cmd.getOptionValue(TUMOR_LOCATION_CSV);

        // Params specific for specific sample
        String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        String purplePurityTsv = cmd.getOptionValue(PURPLE_PURITY_TSV);
        String purpleGeneCnvTsv = cmd.getOptionValue(PURPLE_GENE_CNV_TSV);
        String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);

        if (Utils.anyNull(sampleId,
                knowledgebaseDirectory,
                tumorLocationCsv,
                somaticVariantVcf,
                purplePurityTsv,
                purpleGeneCnvTsv,
                linxFusionTsv,
                cmd.getOptionValue(DB_USER),
                cmd.getOptionValue(DB_PASS),
                cmd.getOptionValue(DB_URL))) {
            printUsageAndExit(options);
        }

        LOGGER.info("Connecting with database");
        DatabaseAccess dbAccess = databaseAccess(cmd);

        LOGGER.info("Reading knowledgebase from {}", knowledgebaseDirectory);
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);

        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationCsv, sampleId);
        List<? extends Variant> passSomaticVariants = readPassSomaticVariants(sampleId, somaticVariantVcf);
        double ploidy = extractPloidy(purplePurityTsv);
        List<GeneCopyNumber> geneCopyNumbers = readGeneCopyNumbers(purpleGeneCnvTsv);
        List<ReportableGeneFusion> geneFusions = readGeneFusions(linxFusionTsv);

        List<EvidenceItem> combinedEvidence = createEvidenceForAllFindings(actionabilityAnalyzer,
                patientPrimaryTumorLocation,
                passSomaticVariants,
                ploidy,
                geneCopyNumbers,
                geneFusions);

        LOGGER.info("Writing evidence items into db");
        dbAccess.writeClinicalEvidence(sampleId, combinedEvidence);
        LOGGER.info("Finished");
    }

    private static double extractPloidy(@NotNull String purplePurityTsv) throws IOException {
        LOGGER.info("Reading purple purity from {}", purplePurityTsv);
        PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info(" Sample ploidy: {}", ploidy);
        return ploidy;
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

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("patient-db - load evidence data", options);
        System.exit(1);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Tumor sample of run");

        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the folder containing knowledgebase files.");
        options.addOption(TUMOR_LOCATION_CSV, true, "Path towards the (curated) tumor location CSV.");

        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_GENE_CNV_TSV, true, "Path towards the purple gene copy number TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");

        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
