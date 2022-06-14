package com.hartwig.hmftools.rose;

import java.io.IOException;

import com.hartwig.hmftools.rose.conclusion.ActionabilityConclusion;
import com.hartwig.hmftools.rose.conclusion.ConclusionAlgo;
import com.hartwig.hmftools.rose.conclusion.RoseConclusionFile;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RoseApplication {

    private static final Logger LOGGER = LogManager.getLogger(RoseApplication.class);
    private static final String VERSION = RoseApplication.class.getPackage().getImplementationVersion();

    public static void main(@NotNull String[] args) throws IOException {
        LOGGER.info("Running ROSE v{}", VERSION);

        Options options = RoseConfig.createOptions();

        RoseConfig config = null;
        try {
            config = RoseConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("ROSE", options);
            System.exit(1);
        }

        new RoseApplication(config).run();

        LOGGER.info("Complete");
    }

    @NotNull
    private final RoseConfig config;

    private RoseApplication(@NotNull final RoseConfig config) {
        this.config = config;
    }

    public void run() throws IOException {
        LOGGER.info("Running ROSE algo on sample {})", config.tumorSampleId());

        RoseAlgo algo = RoseAlgo.build(config.actionabilityDatabaseTsv(),
                config.driverGeneTsv(), config.refGenomeVersion(),
                config.primaryTumorTsv());
        RoseData roseData = algo.run(config);

        ActionabilityConclusion actionabilityConclusion = ConclusionAlgo.generateConclusion(roseData);

        String filename = RoseConclusionFile.generateFilename(config.outputDir(), config.tumorSampleId());
        LOGGER.info("Writing actionability conclusion to file: {}", filename);
        RoseConclusionFile.write(filename, actionabilityConclusion, config.tumorSampleId());
    }
}