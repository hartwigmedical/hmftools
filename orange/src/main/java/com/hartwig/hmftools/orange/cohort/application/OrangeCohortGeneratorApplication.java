package com.hartwig.hmftools.orange.cohort.application;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class OrangeCohortGeneratorApplication {

    private static final Logger LOGGER = LogManager.getLogger(OrangeCohortGeneratorApplication.class);

    private static final String APPLICATION = "ORANGE Cohort Generator";
    public static final String VERSION = OrangeCohortGeneratorApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException, ParseException {
        LOGGER.info("Running {} v{}", APPLICATION, VERSION);

        Options options = OrangeCohortGeneratorConfig.createOptions();

        OrangeCohortGeneratorConfig config = null;
        CommandLine cmd = new DefaultParser().parse(options, args);
        try {
            config = OrangeCohortGeneratorConfig.createConfig(cmd);
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp(APPLICATION, options);
            System.exit(1);
        }

        new OrangeCohortGeneratorApplication(config, cmd).run();
    }

    @NotNull
    private final OrangeCohortGeneratorConfig config;
    @NotNull
    private final CommandLine commandLine;

    public OrangeCohortGeneratorApplication(@NotNull final OrangeCohortGeneratorConfig config, @NotNull final CommandLine commandLine) {
        this.config = config;
        this.commandLine = commandLine;
    }

    private void run() throws IOException {
        LOGGER.info("Generating ORANGE cohorts");

        DatabaseAccess database = DatabaseAccess.createDatabaseAccess(commandLine);

        LOGGER.info("Loading samples from database");
        List<SampleData> samples = SampleDataQuery.run(database);
        LOGGER.info(" Loaded {} samples from database", samples.size());

        LOGGER.info("Loading cohort mappings from {}", config.cohortMappingTsv());
        List<CohortMapping> mappings = CohortMappingFile.read(config.cohortMappingTsv());
        LOGGER.info(" Loaded {} mappings", mappings.size());

        LOGGER.info("Reading DOID model from {}", config.doidJson());
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJson());
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());
        LOGGER.info(" Created DOID parent model from {} edges", doidEntry.edges().size());

        LOGGER.info("Resolving cancer types for {} samples", samples.size());
        CohortMapper mapper = new CohortMapper(doidParentModel, mappings);
        int failed = 0;
        for (SampleData sample : samples) {
            String cancerType = mapper.cancerTypeForDoids(sample.sampleId(), sample.doids());
            if (cancerType.equals(CohortConstants.COHORT_OTHER)) {
                failed++;
            }
        }
        LOGGER.info(" Resolving failed for {} samples", failed);

        LOGGER.info("Done!");
    }
}
