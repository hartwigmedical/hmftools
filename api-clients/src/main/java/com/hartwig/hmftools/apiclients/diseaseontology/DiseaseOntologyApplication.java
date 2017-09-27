package com.hartwig.hmftools.apiclients.diseaseontology;

import com.hartwig.hmftools.apiclients.diseaseontology.api.DiseaseOntologyApiWrapper;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DiseaseOntologyApplication {

    private static final Logger LOGGER = LogManager.getLogger(DiseaseOntologyApplication.class);

    public static void main(final String... args) throws Exception {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final DiseaseOntologyApiWrapper diseaseOntoloy = new DiseaseOntologyApiWrapper();
        diseaseOntoloy.getAllChildrenDoids("DOID:" + 1909).blockingSubscribe(LOGGER::info, LOGGER::error, () -> LOGGER.info("completed"));
        diseaseOntoloy.releaseResources();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
