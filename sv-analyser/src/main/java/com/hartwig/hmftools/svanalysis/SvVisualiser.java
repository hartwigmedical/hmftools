package com.hartwig.hmftools.svanalysis;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.svanalysis.visualisation.CopyNumberAlterations.copyNumberInTracks;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.svanalysis.visualisation.CircosConfigWriter;
import com.hartwig.hmftools.svanalysis.visualisation.CircosDataWriter;
import com.hartwig.hmftools.svanalysis.visualisation.CopyNumberAlteration;
import com.hartwig.hmftools.svanalysis.visualisation.Link;
import com.hartwig.hmftools.svanalysis.visualisation.Segment;
import com.hartwig.hmftools.svanalysis.visualisation.Segments;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SvVisualiser implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(SvVisualiser.class);

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        final Options options = SvVisualiserConfig.createOptions();
        try (final SvVisualiser application = new SvVisualiser(options, args)) {
            application.run();
        } catch (ParseException e) {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SvVisualiser", options);
            System.exit(1);
        }
    }

    private final SvVisualiserConfig config;
    private final ExecutorService executorService;

    private SvVisualiser(final Options options, final String... args) throws ParseException, IOException {
        final CommandLine cmd = createCommandLine(args, options);
        config = SvVisualiserConfig.createConfig(cmd);
        executorService = Executors.newFixedThreadPool(config.threads());
    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        LOGGER.info("Loading data");

        final List<Future<Object>> futures = Lists.newArrayList();
        final List<Integer> clusterIds = config.links().stream().map(Link::clusterId).distinct().sorted().collect(toList());
        for (Integer clusterId : clusterIds) {
            futures.add(executorService.submit(() -> run(clusterId)));
        }

        for (Future<Object> future : futures) {
            future.get();
        }
    }

    @Nullable
    private Object run(int clusterId) throws IOException, InterruptedException {
        final String sample = config.sample() + "." + clusterId;

        final List<Link> clusterLinks = config.links().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Segment> clusterSegments = config.tracks().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Segment> segments = Segments.addMissingTracks(1000, clusterSegments, clusterLinks);
        final List<CopyNumberAlteration> alterations = copyNumberInTracks(100, config.copyNumberAlterations(), segments);

        int maxTracks = segments.stream().mapToInt(Segment::track).max().orElse(0) + 1;
        double maxCopyNumber = alterations.stream().mapToDouble(CopyNumberAlteration::copyNumber).max().orElse(0);
        double maxMinorAllelePloidy = alterations.stream().mapToDouble(CopyNumberAlteration::minorAllelePloidy).max().orElse(0);

        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, config.outputConfPath());
        confWrite.writeConfig(maxTracks, maxCopyNumber, maxMinorAllelePloidy);
        new CircosDataWriter(sample, config.outputConfPath(), maxTracks).write(segments, clusterLinks, alterations);

        final String outputPlotName = sample + ".cluster.png";
        return new CircosExecution(config.circosBin()).generateCircos(confWrite.configPath(), config.outputPlotPath(), outputPlotName);
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close() {
        executorService.shutdown();
        LOGGER.info("Complete");
    }
}
