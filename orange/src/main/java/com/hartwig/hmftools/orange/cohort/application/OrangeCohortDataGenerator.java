package com.hartwig.hmftools.orange.cohort.application;

import static com.hartwig.hmftools.orange.OrangeApplication.APP_NAME;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.cohort.application.OrangeCohortDataGeneratorConfig.createConfig;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.orange.PercentileType;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.orange.cohort.mapping.DoidCohortMapper;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentiles;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentilesFile;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileGenerator;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jetbrains.annotations.NotNull;

public class OrangeCohortDataGenerator
{
    public static void main(String[] args) throws IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        OrangeCohortDataGeneratorConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        OrangeCohortDataGeneratorConfig config = createConfig(configBuilder);

        new OrangeCohortDataGenerator(config, configBuilder).run();
    }

    private final OrangeCohortDataGeneratorConfig config;
    private final DatabaseAccess database;

    public OrangeCohortDataGenerator(final OrangeCohortDataGeneratorConfig config, final ConfigBuilder configBuilder)
    {
        this.config = config;
        database = DatabaseAccess.createDatabaseAccess(configBuilder);
    }

    private void run() throws IOException
    {
        LOGGER.info("Generating ORANGE cohort data");

        LOGGER.info("Loading samples from database");
        List<Sample> samples = SampleQuery.selectFromDatarequest(database);
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
    private static CohortMapper createCohortMapper(@NotNull OrangeCohortDataGeneratorConfig config) throws IOException
    {
        LOGGER.info("Loading cohort mappings from {}", config.cohortMappingTsv());
        List<CohortMapping> mappings = CohortMappingFile.read(config.cohortMappingTsv());
        LOGGER.info(" Loaded {} mappings", mappings.size());

        LOGGER.info("Reading DOID model from {}", config.doidJson());
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(config.doidJson());
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        return new DoidCohortMapper(doidParentModel, mappings);
    }
}
