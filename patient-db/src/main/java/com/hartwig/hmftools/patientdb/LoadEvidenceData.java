package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.ImmutableClinicalTrial;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadEvidenceData {


    private static final Logger LOGGER = LogManager.getLogger(LoadEvidenceData.class);
    private static final String KNOWLEDGEBASE_PATH = "knowledgebase_path";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String RUN_DIR = "run_dir";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final String runDirectoryPath = cmd.getOptionValue(RUN_DIR);
        final String knowledgebase_path = cmd.getOptionValue(KNOWLEDGEBASE_PATH);

        final DatabaseAccess dbAccess = databaseAccess(cmd);

        if (Utils.anyNull(userName, password, databaseUrl, runDirectoryPath)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("patient-db - load evidence data", options);
        } else {
            final File runDirectory = new File(runDirectoryPath);
            if (runDirectory.isDirectory()) {
                final String jdbcUrl = "jdbc:" + databaseUrl;
                final DatabaseAccess dbWriter = new DatabaseAccess(userName, password, jdbcUrl);

                RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory.toPath().toString());
                final String sample = runContext.tumorSample();
                ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebase_path);

                LOGGER.info("sample: " + sample);

                LOGGER.info("Reading clinical Data from DB");
                final String primaryTumorLocation = dbAccess.readTumorLocation(sample);

                LOGGER.info(primaryTumorLocation);

                LOGGER.info("Reading somatic variants from DB");
                final List<EnrichedSomaticVariant> variants = dbAccess.readSomaticVariants(sample);
                LOGGER.info(variants);

                Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant =
                        actionabilityAnalyzer.evidenceForSomaticVariants(variants, primaryTumorLocation);

                final List<EvidenceItem> AllEvidenceItemsVariant = extractAllEvidenceItems(evidencePerVariant);
                final List<ClinicalTrial> allClinicalTrialsVariants = extractAllTrials(AllEvidenceItemsVariant);

                LOGGER.info("Writing evidence items of somatic variants to DB");
                LOGGER.info("Counts of all evidence items: " + AllEvidenceItemsVariant.size());
                dbWriter.writeClinicalEvidence(sample, AllEvidenceItemsVariant);
                LOGGER.info("Writing clinical trials of somatic variants to DB");
                LOGGER.info("Counts of all clinical trials: " + allClinicalTrialsVariants.size());
                dbWriter.writeClinicalTrial(sample, allClinicalTrialsVariants);

                LOGGER.info("Reading gene copy from DB");
                final List<GeneCopyNumber> geneCopyNumber = dbAccess.readGeneCopynumbers(sample);
                LOGGER.info(geneCopyNumber);

                Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                        actionabilityAnalyzer.evidenceForCopyNumbers(geneCopyNumber, primaryTumorLocation);

                final List<EvidenceItem> AllEvidenceItemsGeneCopyNumber = extractAllEvidenceItems(evidencePerGeneCopyNumber);
                final List<ClinicalTrial> allClinicalTrialsGeneCopyNumber = extractAllTrials(AllEvidenceItemsGeneCopyNumber);

                LOGGER.info("Writing evidence items of gene copy numbers to DB");
                LOGGER.info("Counts of all evidence items: " + AllEvidenceItemsGeneCopyNumber.size());
                dbWriter.writeClinicalEvidence(sample, AllEvidenceItemsGeneCopyNumber);
                LOGGER.info("Writing clinical trials of gene copy numbers to DB");
                LOGGER.info("Counts of all clinical trials: " + allClinicalTrialsGeneCopyNumber.size());
                dbWriter.writeClinicalTrial(sample, allClinicalTrialsGeneCopyNumber);

                LOGGER.info("Reading gene fusions from DB");
//                final List<GeneFusion> geneFusions = dbAccess.readGeneFusions(sample);
//
//                Map<GeneFusion, List<EvidenceItem>> evidencePerFusion =
//                        actionabilityAnalyzer.evidenceForFusions(geneFusions, primaryTumorLocation);
//
//                final List<EvidenceItem> AllEvidenceItemsGeneFusions = extractAllEvidenceItems(evidencePerFusion);
//                final List<ClinicalTrial> allClinicalTrialsGeneFusions = extractAllTrials(AllEvidenceItemsGeneFusions);
//
//                LOGGER.info("Writing evidence items of gene fusions to DB");
//                LOGGER.info("Counts of all evidence items: " + AllEvidenceItemsGeneFusions.size());
//                dbWriter.writeClinicalEvidence(sample, AllEvidenceItemsGeneFusions);
//                LOGGER.info("Writing clinical trials of gene fusions to DB");
//                LOGGER.info("Counts of all clinical trials: " + allClinicalTrialsGeneFusions.size());
//                dbWriter.writeClinicalTrial(sample, allClinicalTrialsGeneFusions);

            }
        }
    }

    @NotNull
    public static List<ClinicalTrial> extractAllTrials(@NotNull List<EvidenceItem> evidenceItems) {
        List<ClinicalTrial> trials = Lists.newArrayList();
        for (EvidenceItem evidence : evidenceItems) {
            if (evidence.source().isTrialSource()) {
                trials.add(toClinicalTrial(evidence));
            }
        }
        return trials;
    }

    @NotNull
    private static ClinicalTrial toClinicalTrial(@NotNull EvidenceItem evidenceItem) {
        return ImmutableClinicalTrial.builder()
                .event(evidenceItem.event())
                .acronym(evidenceItem.drug())
                .source(evidenceItem.source())
                .reference(evidenceItem.reference())
                .isOnLabel(evidenceItem.isOnLabel())
                .cancerType(evidenceItem.cancerType())
                .scope(evidenceItem.scope())
                .build();
    }

    @NotNull
    public static List<EvidenceItem> extractAllEvidenceItems(@NotNull Map<?, List<EvidenceItem>> evidenceItemMap) {
        return toList(evidenceItemMap);
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
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(KNOWLEDGEBASE_PATH, true, "Path towards the folder containing knowledgebase files.");
        options.addOption(RUN_DIR, true, "Path towards the folder containing patient run.");
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
