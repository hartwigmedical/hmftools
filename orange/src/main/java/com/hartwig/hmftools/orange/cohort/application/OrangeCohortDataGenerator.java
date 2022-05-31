package com.hartwig.hmftools.orange.cohort.application;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.orange.cohort.mapping.DoidCohortMapper;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentiles;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentilesFile;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileGenerator;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class OrangeCohortDataGenerator {

    private static final Logger LOGGER = LogManager.getLogger(OrangeCohortDataGenerator.class);

    private static final String APPLICATION = "ORANGE Cohort Data Generator";
    private static final String VERSION = OrangeCohortDataGenerator.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws IOException, ParseException {
        LOGGER.info("Running {} v{}", APPLICATION, VERSION);

        Options options = OrangeCohortDataGeneratorConfig.createOptions();

        OrangeCohortDataGeneratorConfig config = null;
        CommandLine cmd = new DefaultParser().parse(options, args);
        try {
            config = OrangeCohortDataGeneratorConfig.createConfig(cmd);
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp(APPLICATION, options);
            System.exit(1);
        }

        new OrangeCohortDataGenerator(config, cmd).run();
    }

    @NotNull
    private final OrangeCohortDataGeneratorConfig config;
    @NotNull
    private final CommandLine commandLine;

    public OrangeCohortDataGenerator(@NotNull final OrangeCohortDataGeneratorConfig config, @NotNull final CommandLine commandLine) {
        this.config = config;
        this.commandLine = commandLine;
    }

    private void run() throws IOException {
        LOGGER.info("Generating ORANGE cohort data");

        DatabaseAccess database = DatabaseAccess.createDatabaseAccess(commandLine);

        LOGGER.info("Loading samples from database");
        List<Sample> samples = SampleQuery.run(database);
        LOGGER.info(" Loaded {} samples from database", samples.size());

        LOGGER.info("Loading SV TMB from database");
        List<Observation> svTmbObservations = SvTmbQuery.run(database, samples);
        LOGGER.info(" Loaded {} SV TMBs from database", svTmbObservations.size());
        database.close();

        PercentileGenerator generator = new PercentileGenerator(createCohortMapper(config));
        Multimap<PercentileType, CohortPercentiles> percentileMap = ArrayListMultimap.create();

        LOGGER.info("Building SV TMB percentiles");
        percentileMap.putAll(PercentileType.SV_TMB, generator.run(svTmbObservations));

        String outputTsv = CohortPercentilesFile.generateOutputTsv(config.outputDirectory());
        LOGGER.info("Writing cohort percentiles to {}", outputTsv);
        CohortPercentilesFile.write(outputTsv, percentileMap);

        LOGGER.info("Done!");
    }

    @NotNull
    private static CohortMapper createCohortMapper(@NotNull OrangeCohortDataGeneratorConfig config) throws IOException {
        LOGGER.info("Loading cohort mappings from {}", config.cohortMappingTsv());
        List<CohortMapping> mappings = CohortMappingFile.read(config.cohortMappingTsv());
        LOGGER.info(" Loaded {} mappings", mappings.size());

        LOGGER.info("Reading DOID model from {}", config.doidJson());
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJson());
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        return new DoidCohortMapper(doidParentModel, mappings);
    }
}
