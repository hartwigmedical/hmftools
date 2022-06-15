package com.hartwig.hmftools.protect;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidenceFile;
import com.hartwig.hmftools.protect.algo.ProtectAlgo;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.ActionableEventsLoader;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProtectApplication {

    private static final Logger LOGGER = LogManager.getLogger(ProtectApplication.class);
    private static final String VERSION = ProtectApplication.class.getPackage().getImplementationVersion();

    public static void main(@NotNull String[] args) throws IOException {
        LOGGER.info("Running PROTECT v{}", VERSION);

        Options options = ProtectConfig.createOptions();

        ProtectConfig config = null;
        try {
            config = ProtectConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PROTECT", options);
            System.exit(1);
        }

        new ProtectApplication(config).run();

        LOGGER.info("Complete");
    }

    @NotNull
    private final ProtectConfig config;

    private ProtectApplication(@NotNull final ProtectConfig config) {
        this.config = config;
    }

    public void run() throws IOException {
        LOGGER.info("Running PROTECT algo on sample {} (with reference sample {})", config.tumorSampleId(), config.referenceSampleId());

        LOGGER.info("Loading DOID file from {}", config.doidJsonFile());
        DoidParents doidParentModel = DoidParents.fromEdges(DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJsonFile()).edges());

        Set<String> patientTumorDoids = patientTumorDoids(config, doidParentModel);
        ActionableEvents actionableEvents = ActionableEventsLoader.readFromDir(config.serveActionabilityDir(), config.refGenomeVersion());

        List<DriverGene> driverGenes = readDriverGenesFromFile(config.driverGeneTsv());

        ProtectAlgo algo = ProtectAlgo.build(actionableEvents, patientTumorDoids, driverGenes);
        List<ProtectEvidence> evidences = algo.run(config);

        String filename = ProtectEvidenceFile.generateFilename(config.outputDir(), config.tumorSampleId());
        LOGGER.info("Writing {} evidence items to file: {}", evidences.size(), filename);
        ProtectEvidenceFile.write(filename, evidences);
    }

    @NotNull
    public static List<DriverGene> readDriverGenesFromFile(@NotNull String driverGeneTsv) throws IOException {
        LOGGER.info(" Reading driver genes from {}", driverGeneTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGeneTsv);
        LOGGER.info("  Read {} driver gene entries", driverGenes.size());
        return driverGenes;
    }

    @NotNull
    private static Set<String> patientTumorDoids(@NotNull ProtectConfig config, @NotNull DoidParents doidParentModel) {
        Set<String> result = Sets.newHashSet();

        Set<String> initialDoids = config.primaryTumorDoids();
        if (initialDoids.isEmpty()) {
            LOGGER.warn("No doids provided for {}. Every treatment will be considered off-label.", config.tumorSampleId());
            return Sets.newHashSet();
        }

        LOGGER.info(" Starting doid resolving for patient with initial tumor doids '{}'", initialDoids);
        for (String initialDoid : initialDoids) {
            result.add(initialDoid);
            result.addAll(doidParentModel.parents(initialDoid));
        }

        LOGGER.info(" {} doids which are considered on-label for patient: '{}'", result.size(), result);
        return result;
    }
}