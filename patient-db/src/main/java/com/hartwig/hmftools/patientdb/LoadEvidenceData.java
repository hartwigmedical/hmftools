package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_PASS;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_URL;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.DB_USER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFile;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.copynumber.ExtractReportableGainsAndLosses;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
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
    private static final String TUMOR_LOCATION_TSV = "tumor_location_tsv";

    private static final String SOMATIC_VARIANT_VCF = "somatic_variant_vcf";
    private static final String PURPLE_DRIVER_CATALOG_TSV = "purple_driver_catalog_tsv";
    private static final String LINX_FUSION_TSV = "linx_fusion_tsv";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createOptions();
        CommandLine cmd = createCommandLine(args, options);

        String sampleId = cmd.getOptionValue(SAMPLE);

        // General params needed for every sample
        String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);
        String tumorLocationTsv = cmd.getOptionValue(TUMOR_LOCATION_TSV);

        // Params specific for specific sample
        String somaticVariantVcf = cmd.getOptionValue(SOMATIC_VARIANT_VCF);
        String purpleDriverCatalogTsv = cmd.getOptionValue(PURPLE_DRIVER_CATALOG_TSV);
        String linxFusionTsv = cmd.getOptionValue(LINX_FUSION_TSV);

        if (Utils.anyNull(sampleId,
                knowledgebaseDirectory,
                tumorLocationTsv,
                somaticVariantVcf,
                purpleDriverCatalogTsv,
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

        String patientPrimaryTumorLocation = extractPatientTumorLocation(tumorLocationTsv, sampleId);
        List<SomaticVariant> passSomaticVariants = readPassSomaticVariants(sampleId, somaticVariantVcf);
        List<ReportableGainLoss> reportableGainLosses = getReportableGainsAndLosses(purpleDriverCatalogTsv);
        List<LinxFusion> geneFusions = readGeneFusions(linxFusionTsv);

        List<EvidenceItem> combinedEvidence = createEvidenceForAllFindings(actionabilityAnalyzer,
                patientPrimaryTumorLocation,
                passSomaticVariants,
                reportableGainLosses,
                geneFusions);

        LOGGER.info("Writing {} evidence items into db for {}", combinedEvidence.size(), sampleId);
        dbAccess.writeClinicalEvidence(sampleId, combinedEvidence);
        LOGGER.info("Complete");
    }

    @NotNull
    private static List<SomaticVariant> readPassSomaticVariants(@NotNull String sampleId, @NotNull String somaticVariantVcf)
            throws IOException {
        LOGGER.info("Reading somatic variants from {}", somaticVariantVcf);
        List<SomaticVariant> passSomaticVariants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sampleId, somaticVariantVcf);
        LOGGER.info(" Loaded {} PASS somatic variants", passSomaticVariants.size());

        List<SomaticVariant> reportedVariants = extractReportedVariants(passSomaticVariants);
        LOGGER.info(" {} variants remaining after filtering on reporting", reportedVariants.size());
        return reportedVariants;
    }

    @NotNull
    private static List<SomaticVariant> extractReportedVariants(@NotNull List<SomaticVariant> somaticVariants) {
        List<SomaticVariant> onlyReportedSomaticVariant = Lists.newArrayList();
        for (SomaticVariant variant : somaticVariants) {
            if (variant.reported()) {
                onlyReportedSomaticVariant.add(variant);
            }
        }
        return onlyReportedSomaticVariant;
    }

    @NotNull
    private static List<ReportableGainLoss> getReportableGainsAndLosses(@NotNull String purpleDriverCatalogTsv) throws IOException {
        LOGGER.info("Reading purple driver catalog from {}", purpleDriverCatalogTsv);
        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(purpleDriverCatalogTsv);
        LOGGER.info(" Loaded {} purple driver catalog records", driverCatalog.size());
        List<ReportableGainLoss> gainsLosses = ExtractReportableGainsAndLosses.toReportableGainsAndLosses(driverCatalog);
        LOGGER.info(" Extracted {} gains and losses from drivers", gainsLosses.size());
        return gainsLosses;
    }

    @NotNull
    private static List<LinxFusion> readGeneFusions(@NotNull String linxFusionTsv) throws IOException {
        LOGGER.info("Reading gene fusions from {}", linxFusionTsv);
        final List<LinxFusion> linxFusions = LinxFusion.read(linxFusionTsv);

        List<LinxFusion> linxFusionsReported = linxFusions.stream().filter(x -> x.reported()).collect(Collectors.toList());

        LOGGER.info(" Loaded {} fusions from {}", linxFusionsReported.size(), linxFusionTsv);
        return linxFusionsReported;
    }

    @NotNull
    private static String extractPatientTumorLocation(@NotNull String tumorLocationTsv, @NotNull String sampleId) throws IOException {
        LOGGER.info("Reading primary tumor locations from {}", tumorLocationTsv);
        List<PatientTumorLocation> patientTumorLocations = PatientTumorLocationFile.read(tumorLocationTsv);
        LOGGER.info(" Loaded tumor locations for {} patients", patientTumorLocations.size());

        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocations, sampleId);

        String patientPrimaryTumorLocation = Strings.EMPTY;
        if (patientTumorLocation != null) {
            patientPrimaryTumorLocation = patientTumorLocation.primaryTumorLocation();
        }

        LOGGER.info(" Retrieved tumor location '{}' for sample {}", patientPrimaryTumorLocation, sampleId);

        return patientPrimaryTumorLocation;
    }

    @NotNull
    private static List<EvidenceItem> createEvidenceForAllFindings(@NotNull ActionabilityAnalyzer actionabilityAnalyzer,
            @NotNull String patientPrimaryTumorLocation, @NotNull List<SomaticVariant> variants,
            @NotNull List<ReportableGainLoss> reportableGainLosses, @NotNull List<LinxFusion> geneFusions) {
        LOGGER.info("Extracting all evidence");

        List<EvidenceItem> evidenceForVariants =
                toList(actionabilityAnalyzer.evidenceForAllVariants(variants, patientPrimaryTumorLocation));
        LOGGER.info(" Found {} evidence items for {} somatic variants.", evidenceForVariants.size(), variants.size());

        List<EvidenceItem> evidenceForGeneCopyNumbers =
                toList(actionabilityAnalyzer.evidenceForCopyNumbers(reportableGainLosses, patientPrimaryTumorLocation));
        LOGGER.info(" Found {} evidence items for {} copy numbers.", evidenceForGeneCopyNumbers.size(), reportableGainLosses.size());

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
        options.addOption(TUMOR_LOCATION_TSV, true, "Path towards the (curated) tumor location TSV.");

        options.addOption(SOMATIC_VARIANT_VCF, true, "Path towards the somatic variant VCF.");
        options.addOption(PURPLE_DRIVER_CATALOG_TSV, true, "Path towards the purple driver catalog TSV.");
        options.addOption(LINX_FUSION_TSV, true, "Path towards the linx fusion TSV.");

        addDatabaseCmdLineArgs(options);
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }
}
