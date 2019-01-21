package com.hartwig.hmftools.svanalysis;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.svanalysis.visualisation.CircosConfigWriter;
import com.hartwig.hmftools.svanalysis.visualisation.CircosDataWriter;
import com.hartwig.hmftools.svanalysis.visualisation.Link;
import com.hartwig.hmftools.svanalysis.visualisation.Track;
import com.hartwig.hmftools.svanalysis.visualisation.Tracks;

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

        LOGGER.info("Extrapolating terminal tracks");
        final List<Track> tracks = Tracks.addLinkTerminals(1000, config.tracks(), config.links());
        final List<Link> links =
                config.links().stream().filter(x -> x.startPosition() != -1 && x.endPosition() != -1).collect(Collectors.toList());

        int maxTracks = tracks.stream().mapToInt(Track::track).max().orElse(0) + 1;

        LOGGER.info("Generating CIRCOS config");
        final CircosConfigWriter confWrite = new CircosConfigWriter(config.sample(), config.outputConfPath(), maxTracks);
        confWrite.writeConfig();
        new CircosDataWriter(config.sample(), config.outputConfPath(), maxTracks).write(tracks, links);

        final String outputPlotName = config.sample() + ".cluster.png";
        new CircosExecution(config.circosBin()).generateCircos(confWrite.configPath(), config.outputPlotPath(), outputPlotName);

    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
