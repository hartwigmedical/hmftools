package com.hartwig.hmftools.svanalysis;

import java.io.IOException;

import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.svanalysis.visualisation.CircosConfigWriter;
import com.hartwig.hmftools.svanalysis.visualisation.CircosDataWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SvVisualiser {

    private static final Logger LOGGER = LogManager.getLogger(SvVisualiser.class);

    public static void main(String[] args) throws IOException, InterruptedException {
        final Options options = SvVisualiserConfig.createOptions();
        try {
            final SvVisualiser application = new SvVisualiser(options, args);
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SvVisualiser", options);
            System.exit(1);
        }
    }

    private final SvVisualiserConfig config;

    private SvVisualiser(final Options options, final String... args) throws ParseException, IOException {
        final CommandLine cmd = createCommandLine(args, options);
        config = SvVisualiserConfig.createConfig(cmd);
    }

    private void run() throws IOException, InterruptedException {

        LOGGER.info("Generating CIRCOS config");
        new CircosDataWriter(config.sample(), config.outputConfPath(), config.maxTracks()).write(config.tracks(), config.links());

        final CircosConfigWriter confWrite = new CircosConfigWriter(config.sample(), config.outputConfPath(), config.maxTracks());
        confWrite.writeConfig();

        final String outputPlotName = config.sample() + ".cluster.png";
        new CircosExecution(config.circosBin()).generateCircos(confWrite.configPath(), config.outputPlotPath(), outputPlotName);

    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
