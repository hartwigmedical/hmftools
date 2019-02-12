package com.hartwig.hmftools.svvisualise;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.svvisualise.circos.Span.allPositions;

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
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.svvisualise.circos.CircosConfigWriter;
import com.hartwig.hmftools.svvisualise.circos.CircosDataWriter;
import com.hartwig.hmftools.svvisualise.circos.ColorPicker;
import com.hartwig.hmftools.svvisualise.circos.Span;
import com.hartwig.hmftools.svvisualise.data.CopyNumberAlteration;
import com.hartwig.hmftools.svvisualise.data.CopyNumberAlterations;
import com.hartwig.hmftools.svvisualise.data.Exon;
import com.hartwig.hmftools.svvisualise.data.Link;
import com.hartwig.hmftools.svvisualise.data.Links;
import com.hartwig.hmftools.svvisualise.data.Segment;
import com.hartwig.hmftools.svvisualise.data.Segments;

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

    private void run() throws InterruptedException, ExecutionException {

        final List<Future<Object>> futures = Lists.newArrayList();
        if (config.singleCluster() != null || config.singleChromosome() != null) {

            if (config.singleCluster() != null) {
                futures.add(executorService.submit(() -> runCluster(config.singleCluster(), false)));
            }

            if (config.singleChromosome() != null) {
                futures.add(executorService.submit(() -> runChromosome(config.singleChromosome())));
            }

        } else {
            final List<Integer> clusterIds = config.links().stream().map(Link::clusterId).distinct().sorted().collect(toList());
            for (Integer clusterId : clusterIds) {
                futures.add(executorService.submit(() -> runCluster(clusterId, true)));
            }

            final Set<String> chromosomes = Sets.newHashSet();
            config.links().stream().map(Link::startChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            config.links().stream().map(Link::endChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            for (final String chromosome : chromosomes) {
                futures.add(executorService.submit(() -> runChromosome(chromosome)));
            }
        }

        for (Future<Object> future : futures) {
            future.get();
        }
    }

    @Nullable
    private Object runChromosome(@NotNull final String chromosome) throws IOException, InterruptedException {
        final String sample = config.sample() + ".chr" + chromosome + (config.debug() ? ".debug" : "");

        final List<Integer> clusterIds = config.links()
                .stream()
                .filter(x -> x.startChromosome().equals(chromosome) || x.endChromosome().equals(chromosome))
                .map(Link::clusterId)
                .collect(toList());

        final List<Link> chromosomeLinks = config.links().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        if (chromosomeLinks.isEmpty()) {
            LOGGER.warn("Chromosome {} not present in file", chromosome);
            return null;
        }

        final List<Segment> chromosomeSegments =
                config.segments().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());

        final List<Exon> chromosomeExons = config.exons().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());

        return runFiltered(sample, chromosomeLinks, chromosomeSegments, chromosomeExons);
    }

    @Nullable
    private Object runCluster(int clusterId, boolean skipSingles) throws IOException, InterruptedException {
        final List<Link> clusterLinks = config.links().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        if (clusterLinks.isEmpty()) {
            LOGGER.warn("Cluster {} not present in file", clusterId);
            return null;
        }

        if (clusterLinks.size() == 1 && skipSingles) {
            LOGGER.info("Skipping simple cluster {}", clusterId);
            return null;
        }

        final String resolvedType = clusterLinks.stream().findFirst().map(Link::resolvedType).orElse("Unknown");

        final String sample =
                config.sample() + ".cluster" + String.format("%03d", clusterId) + "." + resolvedType + ".sv" + clusterLinks.size() + (config
                        .debug() ? ".debug" : "");

        final List<Segment> clusterSegments = config.segments().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Exon> clusterExons = config.exons().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        return runFiltered(sample, clusterLinks, clusterSegments, clusterExons);
    }

    private Object runFiltered(@NotNull final String sample, @NotNull final List<Link> links, @NotNull final List<Segment> filteredSegments,
            @NotNull final List<Exon> exons) throws IOException, InterruptedException {

        final List<GenomePosition> positions = Lists.newArrayList();
        positions.addAll(Links.allPositions(links));
        positions.addAll(allPositions(exons));

        final List<Segment> segments = Segments.ensureCoverage(1000, filteredSegments, links, exons);
        positions.addAll(allPositions(segments));

        final List<CopyNumberAlteration> alterations = CopyNumberAlterations.copyNumbers(100, config.copyNumberAlterations(), Span.span(positions));

        final ColorPicker color = new ColorPicker(links);

        final int chromosomeCount = (int) segments.stream().map(GenomeRegion::chromosome).distinct().count();
        int maxTracks = segments.stream().mapToInt(Segment::track).max().orElse(0) + 1;
        double maxCopyNumber = alterations.stream().mapToDouble(CopyNumberAlteration::copyNumber).max().orElse(0);
        double maxMinorAllelePloidy = alterations.stream().mapToDouble(CopyNumberAlteration::minorAllelePloidy).max().orElse(0);

        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, config.outputConfPath());
        confWrite.writeConfig(chromosomeCount, maxTracks, maxCopyNumber, maxMinorAllelePloidy);
        new CircosDataWriter(config.debug(), color, sample, config.outputConfPath(), maxTracks).write(segments, links, alterations, exons);

        final String outputPlotName = sample + ".png";
        return new CircosExecution(config.circosBin()).generateCircos(confWrite.configPath(),
                config.outputPlotPath(),
                outputPlotName,
                config.outputConfPath());
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
