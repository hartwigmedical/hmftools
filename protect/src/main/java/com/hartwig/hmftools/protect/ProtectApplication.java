package com.hartwig.hmftools.protect;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.bachelor.BachelorDataLoader;
import com.hartwig.hmftools.protect.chord.ChordDataLoader;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
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

        try (ProtectApplication application = new ProtectApplication(options, args)) {
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
        CommandLine cmd = new DefaultParser().parse(options, args);
        this.dbAccess = DatabaseAccess.databaseAccess(cmd);
        this.protectConfig = ProtectConfig.createConfig(cmd);
    }

    public void run() throws IOException {
        List<ProtectEvidence> evidence = protectEvidence(protectConfig);

        LOGGER.info("Writing {} records to database", evidence.size());
        dbAccess.writeProtectEvidence(protectConfig.tumorSampleId(), evidence);

        String filename = ProtectEvidenceFile.generateFilename(protectConfig.outputDir(), protectConfig.tumorSampleId());
        LOGGER.info("Writing {} records to file: {}", evidence.size(), filename);
        ProtectEvidenceFile.write(filename, evidence);
    }

    @NotNull
    private static Set<String> patientTumorDoids(@NotNull ProtectConfig config) throws IOException {
        Set<String> result = Sets.newHashSet();
        LOGGER.info("Loading DOID file from {}", config.doidJsonFile());
        DoidParents doidParentModel = new DoidParents(DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJsonFile()).edges());

        Set<String> initialDoids = config.primaryTumorDoids();
        LOGGER.info("Starting with initial primary tumor doids '{}'", initialDoids);
        for (String initialDoid : initialDoids) {
            result.add(initialDoid);
            result.addAll(doidParentModel.parents(initialDoid));
        }

        LOGGER.info(" Resolved doid tree: {}", result);
        return result;
    }

    @NotNull
    private static List<ProtectEvidence> protectEvidence(@NotNull ProtectConfig config) throws IOException {
        PurpleData purpleData = PurpleDataLoader.load(config);
        LinxData linxData = LinxDataLoader.load(config);
        BachelorData bachelorData = BachelorDataLoader.load(config, purpleData, linxData);
        ChordAnalysis chordAnalysis = ChordDataLoader.load(config);
        Set<String> patientTumorDoids = patientTumorDoids(config);

        ActionableEvents actionableEvents = ActionableEventsLoader.readFromDir(config.serveActionabilityDir(), REF_GENOME_VERSION);

        ProtectAlgo algo = ProtectAlgo.buildAlgoFromServeActionability(actionableEvents, patientTumorDoids);

        return algo.determineEvidence(purpleData, linxData, bachelorData, chordAnalysis);
    }

    @Override
    public void close() {
        dbAccess.close();

        LOGGER.info("Complete");
    }
}
