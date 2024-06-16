package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;
import static java.util.stream.Collectors.toSet;

import static com.hartwig.hmftools.linx.LinxApplication.APP_NAME;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.linx.visualiser.circos.ChromosomeRangeExecution;
import com.hartwig.hmftools.linx.visualiser.circos.CircosConfigWriter;
import com.hartwig.hmftools.linx.visualiser.circos.CircosData;
import com.hartwig.hmftools.linx.visualiser.circos.CircosDataWriter;
import com.hartwig.hmftools.linx.visualiser.circos.ColorPicker;
import com.hartwig.hmftools.linx.visualiser.circos.FusionDataWriter;
import com.hartwig.hmftools.linx.visualiser.circos.FusionExecution;
import com.hartwig.hmftools.linx.visualiser.circos.Span;
import com.hartwig.hmftools.linx.visualiser.data.VisCopyNumbers;
import com.hartwig.hmftools.linx.visualiser.data.VisLinks;
import com.hartwig.hmftools.linx.visualiser.data.VisSegments;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SvVisualiser implements AutoCloseable
{
    public static final Logger VIS_LOGGER = LogManager.getLogger(SvVisualiser.class);

    private final VisualiserConfig mConfig;
    private final SampleData mSampleData;
    private final CircosConfig mCircosConfig;
    private final ExecutorService mExecutorService;

    private final List<Callable<Object>> mCallableImages;
    private final List<Callable<Object>> mCallableConfigs;

    private SvVisualiser(final ConfigBuilder configBuilder) throws Exception
    {
        VIS_LOGGER.info("loading visualiser data");

        mCircosConfig = new CircosConfig(configBuilder);
        mConfig = new VisualiserConfig(configBuilder);
        mSampleData = new SampleData(mConfig);
        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads);

        mCallableImages = Lists.newArrayList();
        mCallableConfigs = Lists.newArrayList();
    }

    private void run() throws Exception
    {
        if(!mCircosConfig.isValid())
        {
            System.exit(1);
        }

        final List<Future<Object>> futures = Lists.newArrayList();

        if(mConfig.PlotReportableEvents)
        {
            Set<Integer> reportableClusterIds = mSampleData.findReportableClusters();

            for(Integer clusterId : reportableClusterIds)
            {
                submitCluster(Lists.newArrayList(clusterId), Collections.EMPTY_LIST, true);
            }
        }
        else if(!mConfig.ClusterIds.isEmpty() || !mConfig.Chromosomes.isEmpty())
        {
            if(!mConfig.ClusterIds.isEmpty())
            {
                submitCluster(mConfig.ClusterIds, mConfig.ChainIds, false);
            }

            if(!mConfig.Chromosomes.isEmpty())
            {
                submitChromosome(mConfig.Chromosomes);
            }
        }
        else
        {
            final List<Integer> clusterIds = mSampleData.SvData.stream().map(x -> x.ClusterId).distinct().sorted().collect(toList());

            for(Integer clusterId : clusterIds)
            {
                submitCluster(Lists.newArrayList(clusterId), Collections.EMPTY_LIST, true);
            }

            final Set<String> chromosomes = Sets.newHashSet();
            mSampleData.SvData.stream().map(x -> x.ChrStart).filter(HumanChromosome::contains).forEach(chromosomes::add);
            mSampleData.SvData.stream().map(x -> x.ChrEnd).filter(HumanChromosome::contains).forEach(chromosomes::add);
            for(final String chromosome : chromosomes)
            {
                submitChromosome(Lists.newArrayList(chromosome));
            }
        }

        mCallableConfigs.forEach(x -> futures.add(mExecutorService.submit(x)));
        mCallableImages.forEach(x -> futures.add(mExecutorService.submit(x)));

        for(Future<Object> future : futures)
        {
            future.get();
        }

        VIS_LOGGER.info("Linx Visualiser complete");
    }

    private void submitChromosome(final List<String> chromosomes)
    {
        if(chromosomes.stream().anyMatch(x -> !HumanChromosome.contains(x)))
        {
            VIS_LOGGER.warn("invalid chromosomes: {}", chromosomes.toString());
            return;
        }

        final String chromosomesStr = chromosomes.size() >= HumanChromosome.values().length
                ? "All"
                : chromosomes.stream().map(RefGenomeFunctions::stripChrPrefix).collect(Collectors.joining("-"));

        final Predicate<VisSvData> linePredicate = x -> !x.isLineElement() || mConfig.IncludeLineElements;
        final Predicate<VisSvData> chromosomePredicate = x -> chromosomes.contains(x.ChrStart) || chromosomes.contains(x.ChrEnd);
        final Predicate<VisSvData> combinedPredicate = chromosomePredicate.and(linePredicate);

        final String sample = mConfig.Sample + ".chr" + chromosomesStr + (mConfig.Debug ? ".debug" : "");

        final Set<Integer> clusterIds = mSampleData.SvData
                .stream()
                .filter(combinedPredicate)
                .map(x -> x.ClusterId)
                .collect(toSet());

        final List<VisSvData> chromosomeLinks = mSampleData.SvData.stream().filter(x -> clusterIds.contains(x.ClusterId)).collect(toList());
        if(chromosomeLinks.isEmpty())
        {
            VIS_LOGGER.warn("chromosomes({}) not present in file", chromosomesStr);
            return;
        }

        if(mCircosConfig.exceedsMaxPlotSvCount(chromosomeLinks.size()))
        {
            VIS_LOGGER.warn("chromosomes({}) svCount({}) exceeds limit({})",
                    chromosomesStr, chromosomeLinks.size(), mCircosConfig.MaxPlotSvCount);
            return;
        }

        final List<VisSegment> chromosomeSegments =
                mSampleData.Segments.stream().filter(x -> clusterIds.contains(x.ClusterId)).collect(toList());

        for(String chromosome : chromosomes)
        {
            chromosomeSegments.add(VisSegments.entireChromosome(mConfig.Sample, chromosome, mConfig.RefGenomeCoords));
        }

        final Set<String> chromosomesOfInterest = Sets.newHashSet(chromosomes);
        chromosomeLinks.forEach(x ->
        {
            chromosomesOfInterest.add(x.ChrStart);
            chromosomesOfInterest.add(x.ChrEnd);
        });
        chromosomeSegments.forEach(x -> chromosomesOfInterest.add(x.chromosome()));

        final List<VisGeneExon> chromosomeExons =
                mSampleData.Exons.stream().filter(x -> chromosomesOfInterest.contains(x.Chromosome)).collect(toList());

        final List<VisProteinDomain> chromosomeProteinDomains =
                mSampleData.ProteinDomains.stream().filter(x -> chromosomesOfInterest.contains(x.chromosome())).collect(toList());

        submitFiltered(ColorPicker::clusterColors, sample, chromosomeLinks, chromosomeSegments, chromosomeExons, chromosomeProteinDomains,
                Collections.emptyList(), false);
    }

    private void submitCluster(final List<Integer> clusterIds, final List<Integer> chainIds, boolean skipSingles)
    {
        final List<VisSvData> clusterSvs = mSampleData.SvData.stream()
                .filter(x -> clusterIds.contains(x.ClusterId))
                .filter(x -> chainIds.isEmpty() || chainIds.contains(x.ChainId))
                .collect(toList());

        String clusterIdsStr = clusterIds.stream().map(x -> String.valueOf(x)).collect(Collectors.joining("-"));

        if(mCircosConfig.exceedsMaxPlotSvCount(clusterSvs.size()))
        {
            VIS_LOGGER.warn("clusterIds({}) svCount({}) exceeds limit({})",
                    clusterIdsStr, clusterSvs.size(), mCircosConfig.MaxPlotSvCount);
            return;
        }

        final List<VisSegment> clusterSegments = mSampleData.Segments.stream()
                .filter(x -> clusterIds.contains(x.ClusterId))
                .filter(x -> chainIds.isEmpty() || chainIds.contains(x.ChainId))
                .collect(toList());

        final List<VisGeneExon> clusterExons =
                mSampleData.Exons.stream().filter(x -> clusterIds.contains(x.ClusterId)).distinct().collect(toList());

        if(clusterSvs.isEmpty())
        {
            VIS_LOGGER.warn("cluster {} not present in file", clusterIdsStr);
            return;
        }

        if(clusterSvs.size() == 1 && skipSingles && clusterExons.isEmpty())
        {
            VIS_LOGGER.debug("skipping simple cluster {}", clusterIdsStr);
            return;
        }

        final Set<Integer> linkChainIds = clusterSvs.stream().map(x -> x.ChainId).collect(Collectors.toSet());
        final Set<Integer> segmentChainIds = clusterSegments.stream().map(x -> x.ChainId).collect(Collectors.toSet());
        segmentChainIds.removeAll(linkChainIds);
        if(!segmentChainIds.isEmpty())
        {
            VIS_LOGGER.warn("Cluster {} contains chain ids {} not found in the links", clusterIdsStr, segmentChainIds);
            return;
        }

        String fileId = mConfig.Sample;

        if(!mConfig.Genes.isEmpty() && mConfig.RestrictClusterByGene)
        {
            StringJoiner genesSj = new StringJoiner("-");
            mConfig.Genes.forEach(x -> genesSj.add(x));

            fileId += mConfig.Genes.size() > 1 ? ".genes-" : ".gene-";
            fileId += genesSj.toString();
        }
        else
        {
            fileId += clusterIds.size() > 1 ? ".clusters-" : ".cluster-";
            fileId += clusterIdsStr;

            if(mConfig.ClusterIds.size() == 1)
            {
                final String resolvedTypeString = clusterSvs.get(0).ClusterResolvedType.toString();
                fileId += "." + resolvedTypeString;
            }
        }

        fileId += ".sv" + clusterSvs.size();

        if(mConfig.Debug)
            fileId += ".debug";

        final List<VisProteinDomain> clusterProteinDomains =
                mSampleData.ProteinDomains.stream().filter(x -> clusterIds.contains(x.ClusterId)).distinct().collect(toList());

        final List<VisFusion> clusterFusions = mSampleData.Fusions.stream().filter(x -> clusterIds.contains(x.ClusterId)).collect(toList());

        submitFiltered(clusterIds.size() == 1 ? ColorPicker::chainColors : ColorPicker::clusterColors,
                fileId, clusterSvs, clusterSegments, clusterExons, clusterProteinDomains, clusterFusions, true);
    }

    private void submitFiltered(final ColorPickerFactory colorPickerFactory,
            final String sample,
            final List<VisSvData> filteredLinks,
            final List<VisSegment> filteredSegments,
            final List<VisGeneExon> filteredExons,
            final List<VisProteinDomain> filteredProteinDomains,
            final List<VisFusion> filteredFusions,
            boolean showSimpleSvSegments)
    {

        final List<GenomePosition> positionsToCover = Lists.newArrayList();
        positionsToCover.addAll(VisLinks.allPositions(filteredLinks));
        positionsToCover.addAll(Span.allPositions(filteredSegments));
        positionsToCover.addAll(Span.allPositions(filteredExons));

        // Limit copy numbers to within segments, links and exons (plus a little extra)
        final List<VisCopyNumber> copyNumbers = VisCopyNumbers.copyNumbers(mSampleData.CopyNumbers, Span.spanPositions(positionsToCover));
        positionsToCover.addAll(Span.allPositions(copyNumbers));

        // Need to extend terminal segments past any current segments, links and exons and copy numbers
        final List<VisSegment> segments = VisSegments.extendTerminals(
                0, filteredSegments, filteredLinks, positionsToCover, showSimpleSvSegments, mConfig.RefGenomeCoords);

        final List<VisSvData> links = VisLinks.addFrame(segments, filteredLinks);

        final ColorPicker color = colorPickerFactory.create(links);

        final CircosData circosData = new CircosData(
                showSimpleSvSegments, mCircosConfig, segments, links, copyNumbers, filteredExons, filteredFusions);

        final CircosConfigWriter confWrite = new CircosConfigWriter(sample, mConfig.OutputConfPath, circosData, mCircosConfig);
        final FusionDataWriter fusionDataWriter = new FusionDataWriter(filteredFusions, filteredExons, filteredProteinDomains);

        mCallableConfigs.add(() -> new CircosDataWriter(color, sample, mConfig.OutputConfPath, mCircosConfig, confWrite, circosData).write());
        if(!fusionDataWriter.finalExons().isEmpty())
        {
            mCallableConfigs.add(() -> fusionDataWriter.write(sample, mConfig.OutputConfPath));
        }

        int minFrame = mCircosConfig.Step ? 0 : circosData.maxFrame();
        for(int frame = minFrame; frame <= circosData.maxFrame(); frame++)
        {
            boolean plotFusion = !fusionDataWriter.finalExons().isEmpty();
            submitFrame(frame, plotFusion, circosData.labelSize(), sample, confWrite);
        }
    }

    private void submitFrame(int frame, boolean hasFusion, double labelSize, String sample, final CircosConfigWriter confWrite)
    {
        boolean plotFusion = hasFusion && !mConfig.Debug;
        boolean plotChromosome = !mConfig.Debug;

        mCallableConfigs.add(() -> confWrite.writeConfig(frame));
        mCallableImages.add(() -> createImageFrame(frame, labelSize, sample, plotFusion, plotChromosome));
    }

    private Object createImageFrame(
            int frame, double labelSize, final String sample, boolean plotFusion, boolean plotChromosome) throws Exception
    {
        final String confFileName = sample + ".circos." + String.format("%03d", frame) + ".conf";
        final String outputFileName = sample + "." + String.format("%03d", frame) + ".png";

        double rLabelSize = 1.2 * labelSize;

        final Object circosResult = new CircosExecution(
                mConfig.CircosBin).generateCircos(mConfig.OutputConfPath + File.separator + confFileName,
                mConfig.OutputPlotPath, outputFileName);

        if(plotFusion)
        {
            int execution = new FusionExecution(sample, outputFileName, mConfig.OutputConfPath, mConfig.OutputPlotPath)
                    .executeR(mCircosConfig, rLabelSize);

            if(execution != 0)
                throw new Exception("plotting error");
        }

        if(plotChromosome)
        {
            int execution = new ChromosomeRangeExecution(
                    sample, mConfig.RefGenVersion, outputFileName, mConfig.OutputConfPath, mConfig.OutputPlotPath)
                    .executeR(mCircosConfig, rLabelSize);

            if(execution != 0)
                throw new Exception("plotting error");
        }

        return circosResult;
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        VisualiserConfig.registerConfig(configBuilder);
        CircosConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        try(final SvVisualiser application = new SvVisualiser(configBuilder))
        {
            application.run();
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.exit(1);
        }
    }

    @Override
    public void close()
    {
        mExecutorService.shutdown();
    }

    private interface ColorPickerFactory
    {
        ColorPicker create(final List<VisSvData> links);
    }
}
