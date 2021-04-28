package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.linx.visualiser.circos.ChromosomeRangeExecution;
import com.hartwig.hmftools.linx.visualiser.circos.CircosConfigWriter;
import com.hartwig.hmftools.linx.visualiser.circos.CircosData;
import com.hartwig.hmftools.linx.visualiser.circos.CircosDataWriter;
import com.hartwig.hmftools.linx.visualiser.circos.ColorPicker;
import com.hartwig.hmftools.linx.visualiser.circos.FusionDataWriter;
import com.hartwig.hmftools.linx.visualiser.circos.FusionExecution;
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

public class SvVisualiser implements AutoCloseable
{
    private static final Logger LOGGER = LogManager.getLogger(SvVisualiser.class);

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException
    {
        final Options options = SvVisualiserConfig.createOptions();
        SvCircosConfig.addOptions(options);

        try (final SvVisualiser application = new SvVisualiser(options, args))
        {
            application.run();
        }
        catch (ParseException e)
        {
            LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SvVisualiser", options);
            System.exit(1);
        }
    }

    private final SvVisualiserConfig config;
    private final SvCircosConfig circosConfig;
    private final ExecutorService executorService;

    private final List<Callable<Object>> callableImages;
    private final List<Callable<Object>> callableConfigs;

    private SvVisualiser(final Options options, final String... args) throws ParseException, IOException
    {
        final CommandLine cmd = createCommandLine(args, options);
        LOGGER.info("Loading data");
        circosConfig = SvCircosConfig.createConfig(cmd);
        config = SvVisualiserConfig.createConfig(cmd);
        executorService = Executors.newFixedThreadPool(config.threads());

        callableImages = Lists.newArrayList();
        callableConfigs = Lists.newArrayList();
    }

    private void run() throws InterruptedException, ExecutionException
    {
        final List<Future<Object>> futures = Lists.newArrayList();

        if (!config.clusters().isEmpty() || !config.chromosomes().isEmpty())
        {
            if (!config.clusters().isEmpty())
            {
                submitCluster(config.clusters(), false);
            }

            if(!config.chromosomes().isEmpty())
            {
                submitChromosome(config.chromosomes());
            }
        }
        else
        {
            final List<Integer> clusterIds = config.links().stream().map(Link::clusterId).distinct().sorted().collect(toList());

            for (Integer clusterId : clusterIds)
            {
                submitCluster(Lists.newArrayList(clusterId), true);
            }

            final Set<String> chromosomes = Sets.newHashSet();
            config.links().stream().map(Link::startChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            config.links().stream().map(Link::endChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            for (final String chromosome : chromosomes)
            {
                submitChromosome(Lists.newArrayList(chromosome));
            }
        }

        callableConfigs.forEach(x -> futures.add(executorService.submit(x)));
        callableImages.forEach(x -> futures.add(executorService.submit(x)));

        for (Future<Object> future : futures)
        {
            future.get();
        }
    }

    private void submitChromosome(@NotNull final List<String> chromosomes)
    {
        if (chromosomes.stream().anyMatch(x -> !HumanChromosome.contains(x)))
        {
            LOGGER.warn("Invalid chromosomes: {}", chromosomes.toString());
            return;
        }

        final String chromosomesStr = chromosomes.size() >= 23
                ? "All"
                : chromosomes.stream().map(RefGenomeFunctions::stripChrPrefix).collect(Collectors.joining("-"));

        final Predicate<Link> linePredicate = x -> !x.isLineElement() || config.includeLineElements();
        final Predicate<Link> chromosomePredicate = x -> chromosomes.contains(x.startChromosome()) || chromosomes.contains(x.endChromosome());
        final Predicate<Link> combinedPredicate = chromosomePredicate.and(linePredicate);

        final String sample = config.sample() + ".chr" + chromosomesStr + (config.debug() ? ".debug" : "");

        final Set<Integer> clusterIds = config.links()
                .stream()
                .filter(combinedPredicate)
                .map(Link::clusterId)
                .collect(toSet());

        final List<Link> chromosomeLinks = config.links().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        if (chromosomeLinks.isEmpty())
        {
            LOGGER.warn("Chromosomes {} not present in file", chromosomesStr);
            return;
        }

        final List<Segment> chromosomeSegments =
                config.segments().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());

        for(String chromosome : chromosomes)
        {
            chromosomeSegments.add(Segments.entireChromosome(config.sample(), chromosome));
        }

        final Set<String> chromosomesOfInterest = Sets.newHashSet(chromosomes);
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

        submitFiltered(ColorPicker::clusterColors, sample, chromosomeLinks, chromosomeSegments, chromosomeExons, chromosomeProteinDomains,
                Collections.emptyList(), false);
    }

    private void submitCluster(final List<Integer> clusterIds, boolean skipSingles)
    {
        final List<Link> clusterLinks = config.links().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        final List<Segment> clusterSegments = config.segments().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        final List<Exon> clusterExons =
                config.exons().stream().filter(x -> clusterIds.contains(x.clusterId())).distinct().collect(toList());

        String clusterIdsStr = "";

        for (int clusterId : clusterIds)
        {
            clusterIdsStr = appendStr(clusterIdsStr, String.valueOf(clusterId), '-');
        }

        if (clusterLinks.isEmpty())
        {
            LOGGER.warn("Cluster {} not present in file", clusterIdsStr);
            return;
        }

        if (clusterLinks.size() == 1 && skipSingles && clusterExons.isEmpty())
        {
            LOGGER.debug("Skipping simple cluster {}", clusterIdsStr);
            return;
        }

        final Set<Integer> linkChainIds = clusterLinks.stream().map(Link::chainId).collect(Collectors.toSet());
        final Set<Integer> segmentChainIds = clusterSegments.stream().map(Segment::chainId).collect(Collectors.toSet());
        segmentChainIds.removeAll(linkChainIds);
        if (!segmentChainIds.isEmpty())
        {
            LOGGER.warn("Cluster {} contains chain ids {} not found in the links", clusterIdsStr, segmentChainIds);
            return;
        }

        final String resolvedTypeString = clusterLinks.stream().findFirst().map(Link::resolvedType).map(Enum::toString).orElse("Unknown");

        final String sample = config.sample() + ".cluster" + clusterIdsStr + "." + resolvedTypeString + ".sv" + clusterLinks.size()
                + (config.debug() ? ".debug" : "");



        final List<ProteinDomain> clusterProteinDomains =
                config.proteinDomain().stream().filter(x -> clusterIds.contains(x.clusterId())).distinct().collect(toList());

        final List<Fusion> clusterFusions = config.fusions().stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());

        submitFiltered(clusterIds.size() == 1 ? ColorPicker::chainColors : ColorPicker::clusterColors,
                sample, clusterLinks, clusterSegments, clusterExons, clusterProteinDomains, clusterFusions, true);
    }

    private void submitFiltered(@NotNull final ColorPickerFactory colorPickerFactory,
            @NotNull final String sample,
            @NotNull final List<Link> filteredLinks,
            @NotNull final List<Segment> filteredSegments,
            @NotNull final List<Exon> filteredExons,
            @NotNull final List<ProteinDomain> filteredProteinDomains,
            @NotNull final List<Fusion> filteredFusions,
            boolean showSimpleSvSegments)
    {

        final List<GenomePosition> positionsToCover = Lists.newArrayList();
        positionsToCover.addAll(Links.allPositions(filteredLinks));
        positionsToCover.addAll(Span.allPositions(filteredSegments));
        positionsToCover.addAll(Span.allPositions(filteredExons));

        // Limit copy numbers to within segments, links and exons (plus a little extra)
        final List<CopyNumberAlteration> alterations =
                CopyNumberAlterations.copyNumbers(config.copyNumberAlterations(), Span.spanPositions(positionsToCover));
        positionsToCover.addAll(Span.allPositions(alterations));

        // Need to extend terminal segments past any current segments, links and exons and copy numbers
        final List<Segment> segments = Segments.extendTerminals(0, filteredSegments, filteredLinks, positionsToCover, showSimpleSvSegments);
        final List<Link> links = Links.addFrame(segments, filteredLinks);

        final ColorPicker color = colorPickerFactory.create(links);

        final CircosData circosData =
                new CircosData(showSimpleSvSegments, circosConfig, segments, links, alterations, filteredExons, filteredFusions);
        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, config.outputConfPath(), circosData, circosConfig);
        final FusionDataWriter fusionDataWriter = new FusionDataWriter(filteredFusions, filteredExons, filteredProteinDomains);

        callableConfigs.add(() -> new CircosDataWriter(color, sample, config.outputConfPath(), circosConfig, confWrite, circosData).write());
        if (!fusionDataWriter.finalExons().isEmpty())
        {
            callableConfigs.add(() -> fusionDataWriter.write(sample, config.outputConfPath()));
        }

        int minFrame = circosConfig.step() ? 0 : circosData.maxFrame();
        for (int frame = minFrame; frame <= circosData.maxFrame(); frame++)
        {
            boolean plotFusion = !fusionDataWriter.finalExons().isEmpty();
            submitFrame(frame, plotFusion, circosData.labelSize(), sample, confWrite);
        }
    }

    private void submitFrame(int frame, boolean fusion, double labelSize, String sample, final CircosConfigWriter confWrite)
    {
        boolean plotFusion = !config.debug() && fusion;
        boolean plotChromosome = !config.debug();

        callableConfigs.add(() -> confWrite.writeConfig(frame));
        callableImages.add(() -> createImageFrame(frame, labelSize, sample, plotFusion, plotChromosome));
    }

    private Object createImageFrame(
            int frame,
            double labelSize,
            @NotNull final String sample,
            boolean plotFusion,
            boolean plotChromosome) throws IOException, InterruptedException
    {

        final String confFileName = sample + ".circos." + String.format("%03d", frame) + ".conf";
        final String outputFileName = sample + "." + String.format("%03d", frame) + ".png";

        double rLabelSize = 1.2 * labelSize;

        final Object circosResult =
                new CircosExecution(config.circosBin()).generateCircos(config.outputConfPath() + File.separator + confFileName,
                        config.outputPlotPath(),
                        outputFileName,
                        config.outputConfPath());

        if (plotFusion)
        {
            new FusionExecution(sample, outputFileName, config.outputConfPath(), config.outputPlotPath()).executeR(circosConfig, rLabelSize);
        }

        if (plotChromosome)
        {
            return new ChromosomeRangeExecution(sample, outputFileName, config.outputConfPath(), config.outputPlotPath()).executeR(circosConfig, rLabelSize);
        }

        return circosResult;
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
