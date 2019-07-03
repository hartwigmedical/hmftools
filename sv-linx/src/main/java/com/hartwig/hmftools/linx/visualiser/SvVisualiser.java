package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.linx.visualiser.circos.CircosConfigWriter;
import com.hartwig.hmftools.linx.visualiser.circos.CircosData;
import com.hartwig.hmftools.linx.visualiser.circos.CircosDataWriter;
import com.hartwig.hmftools.linx.visualiser.circos.ColorPicker;
import com.hartwig.hmftools.linx.visualiser.circos.FusionsDataWriter;
import com.hartwig.hmftools.linx.visualiser.circos.Span;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlteration;
import com.hartwig.hmftools.linx.visualiser.data.CopyNumberAlterations;
import com.hartwig.hmftools.linx.visualiser.data.Exon;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.Link;
import com.hartwig.hmftools.linx.visualiser.data.Links;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.data.Segments;

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

public class SvVisualiser implements AutoCloseable
{

    private static final Logger LOGGER = LogManager.getLogger(SvVisualiser.class);

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException
    {
        final Options options = SvVisualiserConfig.createOptions();
        try (final SvVisualiser application = new SvVisualiser(options, args))
        {
            application.run();
        } catch (ParseException e)
        {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SvVisualiser", options);
            System.exit(1);
        }
    }

    private final SvVisualiserConfig config;
    private final ExecutorService executorService;

    private SvVisualiser(final Options options, final String... args) throws ParseException, IOException
    {
        final CommandLine cmd = createCommandLine(args, options);
        LOGGER.info("Loading data");
        config = SvVisualiserConfig.createConfig(cmd);
        executorService = Executors.newFixedThreadPool(config.threads());
    }

    private void run() throws InterruptedException, ExecutionException
    {

        final List<Future<Object>> futures = Lists.newArrayList();
        if (!config.clusters().isEmpty() || !config.chromosomes().isEmpty())
        {
            config.clusters().forEach(clusterId -> futures.add(executorService.submit(() -> runCluster(clusterId, false))));
            config.chromosomes().forEach(chromosome -> futures.add(executorService.submit(() -> runChromosome(chromosome))));
        }
        else
        {
            final List<Integer> clusterIds = config.links().stream().map(Link::clusterId).distinct().sorted().collect(toList());
            for (Integer clusterId : clusterIds)
            {
                futures.add(executorService.submit(() -> runCluster(clusterId, true)));
            }

            final Set<String> chromosomes = Sets.newHashSet();
            config.links().stream().map(Link::startChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            config.links().stream().map(Link::endChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            for (final String chromosome : chromosomes)
            {
                futures.add(executorService.submit(() -> runChromosome(chromosome)));
            }
        }

        for (Future<Object> future : futures)
        {
            future.get();
        }
    }

    @Nullable
    private Object runChromosome(@NotNull final String chromosome) throws IOException, InterruptedException
    {
        if (!HumanChromosome.contains(chromosome))
        {
            LOGGER.warn("Chromosome {} not permitted", chromosome);
            return null;
        }

        final Predicate<Link> linePredicate = x -> !x.isLineElement() || config.includeLineElements();
        final Predicate<Link> chromosomePredicate = x -> x.startChromosome().equals(chromosome) || x.endChromosome().equals(chromosome);
        final Predicate<Link> combinedPredicate = chromosomePredicate.and(linePredicate);

        final String sample = config.sample() + ".chr" + chromosome + (config.debug() ? ".debug" : "");
        final Set<Integer> clusterIds = config.links()
                .stream()
                .filter(combinedPredicate)
                .map(Link::clusterId)
                .collect(toSet());

        final List<Link> chromosomeLinks = config.links().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        if (chromosomeLinks.isEmpty())
        {
            LOGGER.warn("Chromosome {} not present in file", chromosome);
            return null;
        }

        final List<Segment> chromosomeSegments =
                config.segments().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        chromosomeSegments.add(Segments.entireChromosome(config.sample(), chromosome));

        final Set<String> chromosomesOfInterest = Sets.newHashSet(chromosome);
        chromosomeLinks.forEach(x ->
        {
            chromosomesOfInterest.add(x.startChromosome());
            chromosomesOfInterest.add(x.endChromosome());
        });
        chromosomeSegments.forEach(x -> chromosomesOfInterest.add(x.chromosome()));

        final List<Exon> chromosomeExons =
                config.exons().stream().filter(x -> chromosomesOfInterest.contains(x.chromosome())).collect(toList());

        final List<ProteinDomain> chromosomeProteinDomains =
                config.proteinDomain().stream().filter(x -> chromosomesOfInterest.contains(x.chromosome())).collect(toList());

        return runFiltered(ColorPicker::clusterColors, sample, chromosomeLinks, chromosomeSegments, chromosomeExons, chromosomeProteinDomains, Collections
                .emptyList());
    }

    @Nullable
    private Object runCluster(int clusterId, boolean skipSingles) throws IOException, InterruptedException
    {
        final List<Link> clusterLinks = config.links().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Segment> clusterSegments = config.segments().stream().filter(x -> x.clusterId() == clusterId).collect(toList());

        if (clusterLinks.isEmpty())
        {
            LOGGER.warn("Cluster {} not present in file", clusterId);
            return null;
        }

        if (clusterLinks.size() == 1 && skipSingles)
        {
            LOGGER.info("Skipping simple cluster {}", clusterId);
            return null;
        }

        final Set<Integer> linkChainIds = clusterLinks.stream().map(Link::chainId).collect(Collectors.toSet());
        final Set<Integer> segmentChainIds = clusterSegments.stream().map(Segment::chainId).collect(Collectors.toSet());
        segmentChainIds.removeAll(linkChainIds);
        if (!segmentChainIds.isEmpty())
        {
            LOGGER.warn("Cluster {} contains chain ids {} not found in the links", clusterId, segmentChainIds);
            return null;
        }

        final String resolvedTypeString = clusterLinks.stream().findFirst().map(Link::resolvedType).map(Enum::toString).orElse("Unknown");

        final String sample =
                config.sample() + ".cluster" + String.format("%03d", clusterId) + "." + resolvedTypeString + ".sv" + clusterLinks.size() + (config
                        .debug() ? ".debug" : "");

        final List<Exon> clusterExons = config.exons().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<ProteinDomain> clusterProteinDomains =  config.proteinDomain().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        final List<Fusion> clusterFusions =  config.fusions().stream().filter(x -> x.clusterId() == clusterId).collect(toList());
        return runFiltered(ColorPicker::chainColors, sample, clusterLinks, clusterSegments, clusterExons, clusterProteinDomains, clusterFusions);
    }

    private Object runFiltered(@NotNull final ColorPickerFactory colorPickerFactory, @NotNull final String sample,
            @NotNull final List<Link> links,
            @NotNull final List<Segment> filteredSegments,
            @NotNull final List<Exon> filteredExons,
            @NotNull final List<ProteinDomain> filteredProteinDomains,
            @NotNull final List<Fusion> filteredFusions)
            throws IOException, InterruptedException
    {

        final Set<String> fusionGenes = Sets.newHashSet();
        filteredFusions.forEach(x -> {fusionGenes.add(x.geneUp()); fusionGenes.add(x.geneDown());});


        final List<GenomePosition> positionsToCover = Lists.newArrayList();
        positionsToCover.addAll(Links.allPositions(links));
        positionsToCover.addAll(Span.allPositions(filteredSegments));
        positionsToCover.addAll(Span.allPositions(filteredExons));

        // Need to extend terminal segments past any current segments, links and exons
        final List<Segment> segments = Segments.extendTerminals(1000, filteredSegments, links, positionsToCover);
        positionsToCover.addAll(Span.allPositions(segments));

        // Limit copy numbers to within segments, links and exons (plus a little extra)
        final List<CopyNumberAlteration> alterations =
                CopyNumberAlterations.copyNumbers(100, config.copyNumberAlterations(), Span.span(positionsToCover));

        final ColorPicker color = colorPickerFactory.create(links);

        final CircosData circosData = new CircosData(config.scaleExons(), segments, links, alterations, filteredExons, filteredProteinDomains, filteredFusions);
        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, config.outputConfPath(), circosData);
        confWrite.writeConfig();

        new CircosDataWriter(config.debug(), color, sample, config.outputConfPath(), confWrite).write(circosData);
        new FusionsDataWriter(sample, config.outputConfPath()).write(filteredFusions, filteredExons, filteredProteinDomains);

        final String outputPlotName = sample + ".png";
        return new CircosExecution(config.circosBin()).generateCircos(confWrite.configPath(),
                config.outputPlotPath(),
                outputPlotName,
                config.outputConfPath());
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close()
    {
        executorService.shutdown();
        LOGGER.info("Complete");
    }

    private interface ColorPickerFactory
    {
        @NotNull
        ColorPicker create(@NotNull final List<Link> links);
    }

}
