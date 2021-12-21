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
import com.hartwig.hmftools.linx.visualiser.data.VisCopyNumbers;
import com.hartwig.hmftools.linx.visualiser.data.Fusion;
import com.hartwig.hmftools.linx.visualiser.data.VisSvData;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.data.ProteinDomain;
import com.hartwig.hmftools.linx.visualiser.data.Segment;
import com.hartwig.hmftools.linx.visualiser.data.VisSegments;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

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
    public static final Logger VIS_LOGGER = LogManager.getLogger(SvVisualiser.class);

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException
    {
        final Options options = SvVisualiserConfig.createOptions();
        SampleData.addCmdLineOptions(options);
        CircosConfig.addOptions(options);

        try (final SvVisualiser application = new SvVisualiser(options, args))
        {
            application.run();
        }
        catch (ParseException e)
        {
            VIS_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("SvVisualiser", options);
            System.exit(1);
        }
    }

    private final SvVisualiserConfig mConfig;
    private final SampleData mSampleData;
    private final CircosConfig mCircosConfig;
    private final ExecutorService mExecutorService;

    private final List<Callable<Object>> mCallableImages;
    private final List<Callable<Object>> mCallableConfigs;

    private SvVisualiser(final Options options, final String... args) throws ParseException, IOException
    {
        final CommandLine cmd = createCommandLine(args, options);
        VIS_LOGGER.info("Loading data");

        mCircosConfig = new CircosConfig(cmd);
        mConfig = new SvVisualiserConfig(cmd);
        mSampleData = new SampleData(cmd);
        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads);

        mCallableImages = Lists.newArrayList();
        mCallableConfigs = Lists.newArrayList();
    }

    private void run() throws InterruptedException, ExecutionException
    {
        if(!mCircosConfig.isValid())
        {
            System.exit(1);
        }

        final List<Future<Object>> futures = Lists.newArrayList();

        if(mConfig.PlotReportableEvents)
        {
            Set<Integer> reportableClusterIds = mSampleData.findReportableClusters();

            for (Integer clusterId : reportableClusterIds)
            {
                submitCluster(Lists.newArrayList(clusterId), true);
            }
        }

        else if (!mSampleData.Clusters.isEmpty() || !mSampleData.Chromosomes.isEmpty())
        {
            if (!mSampleData.Clusters.isEmpty())
            {
                submitCluster(mSampleData.Clusters, false);
            }

            if(!mSampleData.Chromosomes.isEmpty())
            {
                submitChromosome(mSampleData.Chromosomes);
            }
        }
        else
        {
            final List<Integer> clusterIds = mSampleData.SvData.stream().map(VisSvData::clusterId).distinct().sorted().collect(toList());

            for (Integer clusterId : clusterIds)
            {
                submitCluster(Lists.newArrayList(clusterId), true);
            }

            final Set<String> chromosomes = Sets.newHashSet();
            mSampleData.SvData.stream().map(VisSvData::startChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            mSampleData.SvData.stream().map(VisSvData::endChromosome).filter(HumanChromosome::contains).forEach(chromosomes::add);
            for (final String chromosome : chromosomes)
            {
                submitChromosome(Lists.newArrayList(chromosome));
            }
        }

        mCallableConfigs.forEach(x -> futures.add(mExecutorService.submit(x)));
        mCallableImages.forEach(x -> futures.add(mExecutorService.submit(x)));

        for (Future<Object> future : futures)
        {
            future.get();
        }
    }

    private void submitChromosome(@NotNull final List<String> chromosomes)
    {
        if (chromosomes.stream().anyMatch(x -> !HumanChromosome.contains(x)))
        {
            VIS_LOGGER.warn("Invalid chromosomes: {}", chromosomes.toString());
            return;
        }

        final String chromosomesStr = chromosomes.size() >= 23
                ? "All"
                : chromosomes.stream().map(RefGenomeFunctions::stripChrPrefix).collect(Collectors.joining("-"));

        final Predicate<VisSvData> linePredicate = x -> !x.isLineElement() || mConfig.IncludeLineElements;
        final Predicate<VisSvData> chromosomePredicate = x -> chromosomes.contains(x.startChromosome()) || chromosomes.contains(x.endChromosome());
        final Predicate<VisSvData> combinedPredicate = chromosomePredicate.and(linePredicate);

        final String sample = mSampleData.Sample + ".chr" + chromosomesStr + (mConfig.Debug ? ".debug" : "");

        final Set<Integer> clusterIds = mSampleData.SvData
                .stream()
                .filter(combinedPredicate)
                .map(VisSvData::clusterId)
                .collect(toSet());

        final List<VisSvData> chromosomeLinks = mSampleData.SvData.stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        if (chromosomeLinks.isEmpty())
        {
            VIS_LOGGER.warn("Chromosomes {} not present in file", chromosomesStr);
            return;
        }

        final List<Segment> chromosomeSegments =
                mSampleData.Segments.stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());

        for(String chromosome : chromosomes)
        {
            chromosomeSegments.add(VisSegments.entireChromosome(mSampleData.Sample, chromosome, mConfig.RefGenomeCoords));
        }

        final Set<String> chromosomesOfInterest = Sets.newHashSet(chromosomes);
        chromosomeLinks.forEach(x ->
        {
            chromosomesOfInterest.add(x.startChromosome());
            chromosomesOfInterest.add(x.endChromosome());
        });
        chromosomeSegments.forEach(x -> chromosomesOfInterest.add(x.chromosome()));

        final List<VisGeneExon> chromosomeExons =
                mSampleData.Exons.stream().filter(x -> chromosomesOfInterest.contains(x.Chromosome)).collect(toList());

        final List<ProteinDomain> chromosomeProteinDomains =
                mSampleData.ProteinDomains.stream().filter(x -> chromosomesOfInterest.contains(x.chromosome())).collect(toList());

        submitFiltered(ColorPicker::clusterColors, sample, chromosomeLinks, chromosomeSegments, chromosomeExons, chromosomeProteinDomains,
                Collections.emptyList(), false);
    }

    private void submitCluster(final List<Integer> clusterIds, boolean skipSingles)
    {
        final List<VisSvData> clusterLinks = mSampleData.SvData.stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        final List<Segment> clusterSegments = mSampleData.Segments.stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());
        final List<VisGeneExon> clusterExons =
                mSampleData.Exons.stream().filter(x -> clusterIds.contains(x.ClusterId)).distinct().collect(toList());

        String clusterIdsStr = "";

        for (int clusterId : clusterIds)
        {
            clusterIdsStr = appendStr(clusterIdsStr, String.valueOf(clusterId), '-');
        }

        if (clusterLinks.isEmpty())
        {
            VIS_LOGGER.warn("Cluster {} not present in file", clusterIdsStr);
            return;
        }

        if (clusterLinks.size() == 1 && skipSingles && clusterExons.isEmpty())
        {
            VIS_LOGGER.debug("Skipping simple cluster {}", clusterIdsStr);
            return;
        }

        final Set<Integer> linkChainIds = clusterLinks.stream().map(VisSvData::chainId).collect(Collectors.toSet());
        final Set<Integer> segmentChainIds = clusterSegments.stream().map(Segment::chainId).collect(Collectors.toSet());
        segmentChainIds.removeAll(linkChainIds);
        if (!segmentChainIds.isEmpty())
        {
            VIS_LOGGER.warn("Cluster {} contains chain ids {} not found in the links", clusterIdsStr, segmentChainIds);
            return;
        }

        final String resolvedTypeString = clusterLinks.stream().findFirst().map(VisSvData::resolvedType).map(Enum::toString).orElse("Unknown");

        final String sample = mSampleData.Sample + ".cluster" + clusterIdsStr + "." + resolvedTypeString + ".sv" + clusterLinks.size()
                + (mConfig.Debug ? ".debug" : "");

        final List<ProteinDomain> clusterProteinDomains =
                mSampleData.ProteinDomains.stream().filter(x -> clusterIds.contains(x.clusterId())).distinct().collect(toList());

        final List<Fusion> clusterFusions = mSampleData.Fusions.stream().filter(x -> clusterIds.contains(x.clusterId())).collect(toList());

        submitFiltered(clusterIds.size() == 1 ? ColorPicker::chainColors : ColorPicker::clusterColors,
                sample, clusterLinks, clusterSegments, clusterExons, clusterProteinDomains, clusterFusions, true);
    }

    private void submitFiltered(final ColorPickerFactory colorPickerFactory,
            final String sample,
            final List<VisSvData> filteredLinks,
            final List<Segment> filteredSegments,
            final List<VisGeneExon> filteredExons,
            final List<ProteinDomain> filteredProteinDomains,
            final List<Fusion> filteredFusions,
            boolean showSimpleSvSegments)
    {

        final List<GenomePosition> positionsToCover = Lists.newArrayList();
        positionsToCover.addAll(VisLinks.allPositions(filteredLinks));
        positionsToCover.addAll(Span.allPositions(filteredSegments));
        positionsToCover.addAll(Span.allPositions(filteredExons));

        // Limit copy numbers to within segments, links and exons (plus a little extra)
        final List<CopyNumberAlteration> alterations =
                VisCopyNumbers.copyNumbers(mSampleData.CopyNumberAlterations, Span.spanPositions(positionsToCover));
        positionsToCover.addAll(Span.allPositions(alterations));

        // Need to extend terminal segments past any current segments, links and exons and copy numbers
        final List<Segment> segments = VisSegments.extendTerminals(
                0, filteredSegments, filteredLinks, positionsToCover, showSimpleSvSegments, mConfig.RefGenomeCoords);

        final List<VisSvData> links = VisLinks.addFrame(segments, filteredLinks);

        final ColorPicker color = colorPickerFactory.create(links);

        final CircosData circosData =
                new CircosData(showSimpleSvSegments, mCircosConfig, segments, links, alterations, filteredExons, filteredFusions);
        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, mConfig.OutputConfPath, circosData, mCircosConfig);
        final FusionDataWriter fusionDataWriter = new FusionDataWriter(filteredFusions, filteredExons, filteredProteinDomains);

        mCallableConfigs.add(() -> new CircosDataWriter(color, sample, mConfig.OutputConfPath, mCircosConfig, confWrite, circosData).write());
        if (!fusionDataWriter.finalExons().isEmpty())
        {
            mCallableConfigs.add(() -> fusionDataWriter.write(sample, mConfig.OutputConfPath));
        }

        int minFrame = mCircosConfig.Step ? 0 : circosData.maxFrame();
        for (int frame = minFrame; frame <= circosData.maxFrame(); frame++)
        {
            boolean plotFusion = !fusionDataWriter.finalExons().isEmpty();
            submitFrame(frame, plotFusion, circosData.labelSize(), sample, confWrite);
        }
    }

    private void submitFrame(int frame, boolean fusion, double labelSize, String sample, final CircosConfigWriter confWrite)
    {
        boolean plotFusion = !mConfig.Debug && fusion;
        boolean plotChromosome = !mConfig.Debug;

        mCallableConfigs.add(() -> confWrite.writeConfig(frame));
        mCallableImages.add(() -> createImageFrame(frame, labelSize, sample, plotFusion, plotChromosome));
    }

    private Object createImageFrame(
            int frame,
            double labelSize,
            final String sample,
            boolean plotFusion,
            boolean plotChromosome) throws IOException, InterruptedException
    {
        final String confFileName = sample + ".circos." + String.format("%03d", frame) + ".conf";
        final String outputFileName = sample + "." + String.format("%03d", frame) + ".png";

        double rLabelSize = 1.2 * labelSize;

        final Object circosResult =
                new CircosExecution(mConfig.CircosBin).generateCircos(mConfig.OutputConfPath + File.separator + confFileName,
                        mConfig.OutputPlotPath,
                        outputFileName,
                        mConfig.OutputConfPath);

        if (plotFusion)
        {
            new FusionExecution(sample, outputFileName, mConfig.OutputConfPath, mConfig.OutputPlotPath).executeR(mCircosConfig, rLabelSize);
        }

        if (plotChromosome)
        {
            return new ChromosomeRangeExecution(sample, outputFileName, mConfig.OutputConfPath, mConfig.OutputPlotPath).executeR(mCircosConfig, rLabelSize);
        }

        return circosResult;
    }

    @NotNull
    private static CommandLine createCommandLine(final String[] args, Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @Override
    public void close()
    {
        mExecutorService.shutdown();
        VIS_LOGGER.info("Complete");
    }

    private interface ColorPickerFactory
    {
        @NotNull
        ColorPicker create(@NotNull final List<VisSvData> links);
    }

}
