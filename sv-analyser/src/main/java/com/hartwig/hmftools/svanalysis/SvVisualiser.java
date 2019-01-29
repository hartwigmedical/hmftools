package com.hartwig.hmftools.svanalysis;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.svanalysis.visualisation.CopyNumberAlterations.copyNumberInTracks;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.svanalysis.visualisation.CircosConfigWriter;
import com.hartwig.hmftools.svanalysis.visualisation.CircosDataWriter;
import com.hartwig.hmftools.svanalysis.visualisation.ColorPicker;
import com.hartwig.hmftools.svanalysis.visualisation.ColorPickerCluster;
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

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException, SQLException {
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

    private SvVisualiser(final Options options, final String... args) throws ParseException, IOException, SQLException {
        final CommandLine cmd = createCommandLine(args, options);
        LOGGER.info("Loading data");
        config = SvVisualiserConfig.createConfig(cmd);
        executorService = Executors.newFixedThreadPool(config.threads());
    }

    private void run() throws InterruptedException, ExecutionException, IOException {

        final List<Future<Object>> futures = Lists.newArrayList();
        final List<Integer> clusterIds = config.links().stream().map(Link::clusterId).distinct().sorted().collect(toList());
        for (Integer clusterId : clusterIds) {
            futures.add(executorService.submit(() -> runCluster(clusterId)));
        }

        final Set<String> chromosomes = Sets.newHashSet();
        config.links().stream().map(Link::startChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
        config.links().stream().map(Link::endChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
        for (final String chromosome : chromosomes) {
            futures.add(executorService.submit(() -> runChromsome(chromosome)));
        }

        for (Future<Object> future : futures) {
            future.get();
        }
    }

    @Nullable
    private Object runChromsome(final String chromosome) throws IOException, InterruptedException {
        final String sample = config.sample() + ".chr" + chromosome + (config.debug() ? ".debug" : "");

        final List<Integer> clusterIds = config.links()
                .stream()
                .filter(x -> x.startChromosome().equals(chromosome) || x.endChromosome().equals(chromosome))
                .map(Link::clusterId)
                .collect(toList());

        final List<Link> links = config.links().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        final List<Segment> clusterSegments = config.segments().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        final List<Segment> segments = Segments.extendTerminals(1000, clusterSegments, links);
        final List<CopyNumberAlteration> alterations = copyNumberInTracks(100, config.copyNumberAlterations(), segments);
        final ColorPicker color = new ColorPickerCluster(links);

        final int chromosomeCount = (int) segments.stream().map(GenomeRegion::chromosome).distinct().count();
        final int maxTracks = segments.stream().mapToInt(Segment::track).max().orElse(0) + 1;
        final double maxCopyNumber = alterations.stream().mapToDouble(CopyNumberAlteration::copyNumber).max().orElse(0);
        final double maxMinorAllelePloidy = alterations.stream().mapToDouble(CopyNumberAlteration::minorAllelePloidy).max().orElse(0);

        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, config.outputConfPath());
        confWrite.writeConfig(chromosomeCount, maxTracks, maxCopyNumber, maxMinorAllelePloidy);
        new CircosDataWriter(config.debug(), color, sample, config.outputConfPath(), maxTracks).write(segments, links, alterations);

        final String outputPlotName = sample + ".png";
        return new CircosExecution(config.circosBin()).generateCircos(confWrite.configPath(), config.outputPlotPath(), outputPlotName);
    }

    @Nullable
    private Object runCluster(int clusterId) throws IOException, InterruptedException {
        final String sample = config.sample() + ".cluster" + String.format("%03d", clusterId) + (config.debug() ? ".debug" : "");

        final List<Link> clusterLinks = config.links().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Segment> clusterSegments = config.segments().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Segment> segments = Segments.extendTerminals(1000, clusterSegments, clusterLinks);
        final List<CopyNumberAlteration> alterations = copyNumberInTracks(100, config.copyNumberAlterations(), segments);
        final ColorPicker color = new ColorPicker() {
        };

        final int chromosomeCount = (int) segments.stream().map(GenomeRegion::chromosome).distinct().count();
        int maxTracks = segments.stream().mapToInt(Segment::track).max().orElse(0) + 1;
        double maxCopyNumber = alterations.stream().mapToDouble(CopyNumberAlteration::copyNumber).max().orElse(0);
        double maxMinorAllelePloidy = alterations.stream().mapToDouble(CopyNumberAlteration::minorAllelePloidy).max().orElse(0);

        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, config.outputConfPath());
        confWrite.writeConfig(chromosomeCount, maxTracks, maxCopyNumber, maxMinorAllelePloidy);
        new CircosDataWriter(config.debug(), color, sample, config.outputConfPath(), maxTracks).write(segments, clusterLinks, alterations);

        final String outputPlotName = sample + ".png";
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
