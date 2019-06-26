package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
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
    private static final String KNOWLEDGEBASE_DIRECTORY = "knowledgebase_dir";

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
        final String knowledgebaseDirectory = cmd.getOptionValue(KNOWLEDGEBASE_DIRECTORY);

        if (Utils.anyNull(userName, password, databaseUrl, runDirectoryPath, knowledgebaseDirectory)) {
            printUsageAndExit(options);
        }

        final File runDirectory = new File(runDirectoryPath);
        if (!runDirectory.isDirectory()) {
            LOGGER.warn("run_dir {} has to be an actual directory", runDirectory);
            printUsageAndExit(options);
        }

        DatabaseAccess dbAccess = databaseAccess(cmd);

        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDirectory);

        RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDirectory.toPath().toString());
        String sample = runContext.tumorSample();

        LOGGER.info("sample: " + sample);

        LOGGER.info("Reading primary tumor location from DB");
        String primaryTumorLocation = dbAccess.readTumorLocation(sample);
        LOGGER.info("Primary tumor location: " + primaryTumorLocation);

        LOGGER.info("Reading somatic variants from DB");
        List<EnrichedSomaticVariant> variants = dbAccess.readSomaticVariants(sample);
        List<EnrichedSomaticVariant> passSomaticVariants = extractPassSomaticVariants(variants);

        LOGGER.info("All somatic Variants: " + variants.size());
        LOGGER.info("PASS somatic Variants: " + passSomaticVariants.size());

        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant =
                actionabilityAnalyzer.evidenceForSomaticVariants(passSomaticVariants, primaryTumorLocation);

        List<EvidenceItem> allEvidenceForSomaticVariants = extractAllEvidenceItems(evidencePerVariant);
        LOGGER.info("Found {} evidence items for {} somatic variants.", allEvidenceForSomaticVariants.size(), passSomaticVariants.size());

        LOGGER.info("Reading gene copy numbers and sample ploidy from DB");
        List<GeneCopyNumber> geneCopyNumbers = dbAccess.readGeneCopynumbers(sample);
        LOGGER.info("All geneCopyNumbers: " + geneCopyNumbers.size());

        PurityContext purityContext = dbAccess.readPurityContext(sample);
        assert purityContext != null;

        double ploidy = purityContext.bestFit().ploidy();
        LOGGER.info("Sample ploidy: " + ploidy);

        Map<GeneCopyNumber, List<EvidenceItem>> evidencePerGeneCopyNumber =
                actionabilityAnalyzer.evidenceForCopyNumbers(geneCopyNumbers, primaryTumorLocation, ploidy);

        List<EvidenceItem> allEvidenceForCopyNumbers = extractAllEvidenceItems(evidencePerGeneCopyNumber);
        LOGGER.info("Found {} evidence items for {} copy numbers.",
                allEvidenceForCopyNumbers.size(),
                evidencePerGeneCopyNumber.keySet().size());

        LOGGER.info("Reading gene fusions from DB");
        List<ReportableGeneFusion> fusions = dbAccess.readGeneFusions(sample);
        LOGGER.info(" All fusions: " + fusions.size());

        Map<ReportableGeneFusion, List<EvidenceItem>> evidencePerFusion =
                actionabilityAnalyzer.evidenceForFusions(fusions, primaryTumorLocation);

        List<EvidenceItem> allEvidenceForGeneFusions = extractAllEvidenceItems(evidencePerFusion);
        LOGGER.info("Found {} evidence items for {} gene fusions.", allEvidenceForGeneFusions.size(), fusions.size());

        List<EvidenceItem> combinedEvidence = Lists.newArrayList();
        combinedEvidence.addAll(allEvidenceForSomaticVariants);
        combinedEvidence.addAll(allEvidenceForCopyNumbers);
        combinedEvidence.addAll(allEvidenceForGeneFusions);

        dbAccess.writeClinicalEvidence(sample, combinedEvidence);
    }

    @NotNull
    private static List<EnrichedSomaticVariant> extractPassSomaticVariants(@NotNull List<EnrichedSomaticVariant> somaticVariant) {
        List<EnrichedSomaticVariant> passSomaticVariants = Lists.newArrayList();
        for (EnrichedSomaticVariant variant : somaticVariant) {
            if (!variant.isFiltered()) {
                passSomaticVariants.add(variant);
            }
        }
        return passSomaticVariants;
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
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(KNOWLEDGEBASE_DIRECTORY, true, "Path towards the folder containing knowledgebase files.");
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
