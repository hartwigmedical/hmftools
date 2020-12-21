package com.hartwig.hmftools.protect;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItemFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.bachelor.BachelorDataLoader;
import com.hartwig.hmftools.protect.evidence.CopyNumberEvidence;
import com.hartwig.hmftools.protect.evidence.FusionEvidence;
import com.hartwig.hmftools.protect.evidence.VariantEvidence;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.ActionableEventsLoader;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProtectApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(ProtectApplication.class);

    private static final RefGenomeVersion REF_GENOME_VERSION = RefGenomeVersion.V37;

    public static void main(@NotNull String[] args) throws IOException {
        Options options = ProtectConfig.createOptions();
        DatabaseAccess.addDatabaseCmdLineArgs(options);

        try (final ProtectApplication application = new ProtectApplication(options, args)) {
            application.run();
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PROTECT", options);
            System.exit(1);
        } catch (SQLException exception) {
            LOGGER.warn(exception);
            System.exit(1);
        }
    }

    @NotNull
    private final DatabaseAccess dbAccess;
    @NotNull
    private final ProtectConfig protectConfig;

    public ProtectApplication(final Options options, final String... args) throws ParseException, SQLException, IOException {
        final CommandLine cmd = new DefaultParser().parse(options, args);
        this.dbAccess = DatabaseAccess.databaseAccess(cmd);
        this.protectConfig = ProtectConfig.createConfig(cmd);
    }

    public void run() throws IOException {
        final List<ProtectEvidenceItem> evidence = protectEvidence(protectConfig);

        LOGGER.info("Writing {} records to database", evidence.size());
        dbAccess.writeProtectEvidence(protectConfig.tumorSampleId(), evidence);

        final String filename = ProtectEvidenceItemFile.generateFilename(protectConfig.outputDir(), protectConfig.tumorSampleId());
        LOGGER.info("Writing {} records to file: {}", evidence.size(), filename);
        ProtectEvidenceItemFile.write(filename, evidence);
    }

    @NotNull
    private static Set<String> doids(@NotNull ProtectConfig config) throws IOException {
        final Set<String> result = Sets.newHashSet();
        LOGGER.info("Loading DOID file from {}", config.doidJsonFile());
        final DoidParents doidParent = new DoidParents(DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJsonFile()).edges());

        final Set<String> initialDoids = config.primaryTumorDoids();
        LOGGER.info("Starting with initial primary tumor doids '{}'", initialDoids);
        for (String initialDoid : initialDoids) {
            result.add(initialDoid);
            result.addAll(doidParent.parents(initialDoid));
        }

        LOGGER.info(" Resolved doid tree: {}", String.join(";", result));
        return result;
    }

    @NotNull
    private static List<ProtectEvidenceItem> protectEvidence(@NotNull ProtectConfig config) throws IOException {
        final Set<String> doids = doids(config);

        // Serve Data
        final ActionableEvents actionableEvents = ActionableEventsLoader.readFromDir(config.serveActionabilityDir(), REF_GENOME_VERSION);
        final VariantEvidence variantEvidenceFactory = new VariantEvidence(actionableEvents.hotspots(), actionableEvents.ranges());
        final CopyNumberEvidence copyNumberEvidenceFactory = new CopyNumberEvidence(actionableEvents.genes());
        final FusionEvidence fusionEvidenceFactory = new FusionEvidence(actionableEvents.genes(), actionableEvents.fusions());

        // Additional configuration
        final GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(config.germlineReportingTsv());

        // External Data
        final LinxData linxData = LinxDataLoader.load(config);
        final PurpleData purpleData = PurpleDataLoader.load(config);
        final BachelorData bachelorData = BachelorDataLoader.load(config.bachelorTsv(), purpleData, linxData, germlineReportingModel);

        // Evidence
        final List<ProtectEvidenceItem> variantEvidence =
                variantEvidenceFactory.evidence(doids, bachelorData.germlineVariants(), purpleData.somaticVariants());
        final List<ProtectEvidenceItem> copyNumberEvidence = copyNumberEvidenceFactory.evidence(doids, purpleData.copyNumberAlterations());
        final List<ProtectEvidenceItem> fusionEvidence = fusionEvidenceFactory.evidence(doids, linxData.fusions());

        final List<ProtectEvidenceItem> result = Lists.newArrayList();
        result.addAll(variantEvidence);
        result.addAll(copyNumberEvidence);
        result.addAll(fusionEvidence);
        return result;
    }

    @Override
    public void close() {
        dbAccess.close();
        LOGGER.info("Complete");
    }
}
