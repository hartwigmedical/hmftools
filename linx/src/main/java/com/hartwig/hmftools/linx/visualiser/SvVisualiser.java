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
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.circos.CircosExecution;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.purple.PurpleSegment;
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
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
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
import org.jetbrains.annotations.Nullable;

public class SvVisualiser implements AutoCloseable
{
    public static final Logger VIS_LOGGER = LogManager.getLogger(SvVisualiser.class);

    private final VisualiserConfig mConfig;
    private final SampleData mSampleData;

    @Nullable
    private final EnsemblDataCache mEnsemblDataCache;

    private final CircosConfig mCircosConfig;
    private final ExecutorService mExecutorService;

    private final List<Callable<Object>> mCallableImages;
    private final List<Callable<Object>> mCallableConfigs;

    private SvVisualiser(final ConfigBuilder configBuilder) throws Exception
    {
        mCircosConfig = new CircosConfig(configBuilder);
        mConfig = new VisualiserConfig(configBuilder);
        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads);

        mCallableImages = Lists.newArrayList();
        mCallableConfigs = Lists.newArrayList();

        if(mConfig.EnsemblDataDir != null)
        {
            VIS_LOGGER.info("loading Ensembl data");
            mEnsemblDataCache = new EnsemblDataCache(mConfig.EnsemblDataDir, mConfig.RefGenVersion);
            mEnsemblDataCache.setRequiredData(true, false, false, true);
            mEnsemblDataCache.load(false);
        }
        else
        {
            mEnsemblDataCache = null;
        }

        mSampleData = new SampleData(mConfig);
    }

    private void run() throws Exception
    {
        if(!mCircosConfig.isValid())
        {
            System.exit(1);
        }

        List<Future<Object>> futures = Lists.newArrayList();

        if(mConfig.PlotReportableEvents)
        {
            submitReportableEvents();
        }
        else if(!mConfig.ClusterIds.isEmpty() || !mConfig.Chromosomes.isEmpty() || !mConfig.Genes.isEmpty())
        {
            submitSelectedEvents();
        }
        else
        {
            submitDefaultEvents();
        }

        mCallableConfigs.forEach(x -> futures.add(mExecutorService.submit(x)));
        mCallableImages.forEach(x -> futures.add(mExecutorService.submit(x)));

        for(Future<Object> future : futures)
        {
            future.get();
        }

        VIS_LOGGER.info("Linx Visualiser complete");
    }

    private void submitDefaultEvents()
    {
        List<Integer> clusterIds = mSampleData.SvData.stream().map(x -> x.ClusterId).distinct().sorted().collect(toList());

        for(Integer clusterId : clusterIds)
        {
            submitCluster(Lists.newArrayList(clusterId), Collections.EMPTY_LIST, true);
        }

        Set<String> chromosomes = Sets.newHashSet();
        mSampleData.SvData.stream().map(x -> x.ChrStart).filter(HumanChromosome::contains).forEach(chromosomes::add);
        mSampleData.SvData.stream().map(x -> x.ChrEnd).filter(HumanChromosome::contains).forEach(chromosomes::add);
        for(String chromosome : chromosomes)
        {
            submitChromosome(Lists.newArrayList(chromosome));
        }

        List<DriverCatalog> genesWithCNVs = mSampleData.findGenesWithCNVs();
        if(!genesWithCNVs.isEmpty())
        {
            for(DriverCatalog gene : genesWithCNVs)
            {
                String geneName = gene.gene();
                DriverType geneDriverType = gene.driver();
                String geneChromosome = mEnsemblDataCache.getGeneDataByName(geneName).Chromosome;
                submitChromosome(Lists.newArrayList(geneChromosome), geneName, geneDriverType);
            }
        }

    }

    private void submitReportableEvents()
    {
        Set<Integer> reportableClusterIds = mSampleData.findReportableClusters();

        for(Integer clusterId : reportableClusterIds)
        {
            submitCluster(Lists.newArrayList(clusterId), Collections.EMPTY_LIST, true);
        }
    }

    private void submitSelectedEvents()
    {
        if(!mConfig.ClusterIds.isEmpty())
        {
            submitCluster(mConfig.ClusterIds, mConfig.ChainIds, false);
        }

        if(!mConfig.Chromosomes.isEmpty())
        {
            submitChromosome(mConfig.Chromosomes);
        }

        if(!mConfig.Genes.isEmpty())
        {
            List<DriverCatalog> genesWithCNVs = mSampleData.findGenesWithCNVs();

            for(String geneName : mConfig.Genes)
            {
                String geneChromosome = mEnsemblDataCache.getGeneDataByName(geneName).Chromosome;

                DriverCatalog matchedGeneCnvEntry = genesWithCNVs.stream()
                        .filter(x -> x.gene().equals(geneName))
                        .findFirst()
                        .orElse(null);

                DriverType geneDriverType = null;
                if(matchedGeneCnvEntry != null)
                {
                    geneDriverType = matchedGeneCnvEntry.driver();
                    VIS_LOGGER.debug("found CNV entry({}) in PURPLE driver catalog for manually selected gene({})",
                            geneDriverType, geneName);
                }

                submitChromosome(Lists.newArrayList(geneChromosome), geneName, geneDriverType);
            }
        }
    }

    private void submitChromosome(List<String> chromosomes)
    {
        submitChromosome(chromosomes, null, null);
    }

    private void submitChromosome(List<String> chromosomes, @Nullable final String geneName, @Nullable final DriverType geneDriverType)
    {
        if(chromosomes.stream().anyMatch(x -> !HumanChromosome.contains(x)))
        {
            VIS_LOGGER.warn("invalid chromosomes: {}", chromosomes.toString());
            return;
        }

        String chromosomesStr = chromosomes.size() >= HumanChromosome.values().length
                ? "All"
                : chromosomes.stream().map(RefGenomeFunctions::stripChrPrefix).collect(Collectors.joining("-"));

        StringJoiner fileIdSj = new StringJoiner(".");

        fileIdSj.add(mConfig.Sample);
        fileIdSj.add("chr" + chromosomesStr);

        if(geneName != null)
        {
            fileIdSj.add(geneName);

            if(geneDriverType != null)
                fileIdSj.add(geneDriverType.toString().toLowerCase());
        }

        if(mConfig.Debug)
            fileIdSj.add("debug");

        String fileId = fileIdSj.toString();

        Predicate<VisSvData> linePredicate = x -> !x.isLineElement() || mConfig.IncludeLineElements;
        Predicate<VisSvData> chromosomePredicate = x -> chromosomes.contains(x.ChrStart) || chromosomes.contains(x.ChrEnd);
        Predicate<VisSvData> combinedPredicate = chromosomePredicate.and(linePredicate);

        Set<Integer> clusterIds = mSampleData.SvData
                .stream()
                .filter(combinedPredicate)
                .map(x -> x.ClusterId)
                .collect(toSet());

        List<VisSvData> chromosomeLinks = mSampleData.SvData.stream().filter(x -> clusterIds.contains(x.ClusterId)).collect(toList());

        if(mCircosConfig.exceedsMaxPlotSvCount(chromosomeLinks.size()))
        {
            VIS_LOGGER.warn("plot({}) - chromosomes({}) svCount({}) exceeds limit({})",
                    fileId, chromosomesStr, chromosomeLinks.size(), mCircosConfig.MaxPlotSvCount);
            return;
        }

        List<VisSegment> chromosomeSegments =
                mSampleData.Segments.stream().filter(x -> clusterIds.contains(x.ClusterId)).collect(toList());

        for(String chromosome : chromosomes)
        {
            chromosomeSegments.add(VisSegments.entireChromosome(mConfig.Sample, chromosome, mConfig.RefGenomeCoords));
        }

        Set<String> chromosomesOfInterest = Sets.newHashSet(chromosomes);
        chromosomeLinks.forEach(x ->
        {
            chromosomesOfInterest.add(x.ChrStart);
            chromosomesOfInterest.add(x.ChrEnd);
        });
        chromosomeSegments.forEach(x -> chromosomesOfInterest.add(x.chromosome()));

        List<VisGeneExon> chromosomeExons;
        if(geneName == null)
        {
            chromosomeExons = mSampleData.Exons.stream().filter(x -> chromosomesOfInterest.contains(x.Chromosome)).collect(toList());

            if(!mConfig.Genes.isEmpty())
                chromosomeExons.addAll(getExonDataFromCache(mConfig.Genes, mConfig.ClusterIds, chromosomeExons, fileId));
        }
        else
        {
            chromosomeExons = getExonDataFromCache(Set.of(geneName), mConfig.ClusterIds, Lists.newArrayList(), fileId);
        }

        List<VisProteinDomain> chromosomeProteinDomains =
                mSampleData.ProteinDomains.stream().filter(x -> chromosomesOfInterest.contains(x.chromosome())).collect(toList());

        submitFiltered(ColorPicker::clusterColors, fileId, chromosomeLinks, chromosomeSegments, chromosomeExons, chromosomeProteinDomains,
                Collections.emptyList(), false);
    }

    private List<VisGeneExon> getExonDataFromCache(
            final Set<String> geneList, final List<Integer> clusterIds, final List<VisGeneExon> currentExons, final String fileId
    )
    {
        final List<VisGeneExon> exonList = Lists.newArrayList();

        final List<Integer> allClusterIds = clusterIds.isEmpty() ? Lists.newArrayList(0) : clusterIds;

        for(final String geneName : geneList)
        {
            if(currentExons.stream().anyMatch(x -> x.Gene.equals(geneName) && clusterIds.contains(x.ClusterId)))
                continue;

            VIS_LOGGER.debug("plot({}) - loading exon data for gene({})", fileId, geneName);

            GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);
            TranscriptData transcriptData = geneData != null ? mEnsemblDataCache.getCanonicalTranscriptData(geneData.GeneId) : null;

            if(transcriptData == null)
            {
                VIS_LOGGER.warn("data not found for specified gene({})", geneName);
                continue;
            }

            for(Integer clusterId : allClusterIds)
            {
                exonList.addAll(VisExons.extractExonList(mConfig.Sample, clusterId, geneData, transcriptData));
            }
        }

        return exonList;
    }

    private void submitCluster(final List<Integer> clusterIds, final List<Integer> chainIds, boolean skipSingles)
    {
        List<VisSvData> clusterSvs = mSampleData.SvData.stream()
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

        List<VisSegment> clusterSegments = mSampleData.Segments.stream()
                .filter(x -> clusterIds.contains(x.ClusterId))
                .filter(x -> chainIds.isEmpty() || chainIds.contains(x.ChainId))
                .collect(toList());

        List<VisGeneExon> clusterExons =
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

        Set<Integer> linkChainIds = clusterSvs.stream().map(x -> x.ChainId).collect(Collectors.toSet());
        Set<Integer> segmentChainIds = clusterSegments.stream().map(x -> x.ChainId).collect(Collectors.toSet());
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
                String resolvedTypeString = clusterSvs.get(0).ClusterResolvedType.toString();
                fileId += "." + resolvedTypeString;
            }
        }

        fileId += ".sv" + clusterSvs.size();

        if(mConfig.Debug)
            fileId += ".debug";

        List<VisProteinDomain> clusterProteinDomains =
                mSampleData.ProteinDomains.stream().filter(x -> clusterIds.contains(x.ClusterId)).distinct().collect(toList());

        List<VisFusion> clusterFusions = mSampleData.Fusions.stream().filter(x -> clusterIds.contains(x.ClusterId)).collect(toList());

        submitFiltered(clusterIds.size() == 1 ? ColorPicker::chainColors : ColorPicker::clusterColors,
                fileId, clusterSvs, clusterSegments, clusterExons, clusterProteinDomains, clusterFusions, true);
    }

    private void submitFiltered(final ColorPickerFactory colorPickerFactory,
            final String fileId,
            final List<VisSvData> filteredLinks,
            final List<VisSegment> filteredSegments,
            final List<VisGeneExon> filteredExons,
            final List<VisProteinDomain> filteredProteinDomains,
            final List<VisFusion> filteredFusions,
            boolean showSimpleSvSegments)
    {
        List<GenomePosition> positionsToCover = Lists.newArrayList();
        positionsToCover.addAll(VisLinks.allPositions(filteredLinks));
        positionsToCover.addAll(Span.allPositions(filteredSegments));
        positionsToCover.addAll(Span.allPositions(filteredExons));

        // Limit copy numbers to within segments, links and exons (plus a little extra)
        List<VisCopyNumber> copyNumbers = VisCopyNumbers.copyNumbers(mSampleData.CopyNumbers, Span.spanPositions(positionsToCover));
        positionsToCover.addAll(Span.allPositions(copyNumbers));

        positionsToCover = positionsToCover.stream().distinct().collect(toList());
        Collections.sort(positionsToCover);

        List<GenomeRegion> regionsToCover = Span.spanPositions(positionsToCover);
        List<AmberBAF> filteredAmberBAFs = Lists.newArrayList();
        List<CobaltRatio> filteredCobaltRatios = Lists.newArrayList();
        List<PurpleSegment> filteredPurpleSegments = Lists.newArrayList();
        for(GenomeRegion region : regionsToCover)
        {
            if(mConfig.AmberDir != null)
            {
                List<AmberBAF> regionAmberBAFs = mSampleData.AmberBAFs.stream()
                        .filter(x -> x.Chromosome.equals(region.chromosome()) && region.start()<=x.position() && region.end()>=x.position())
                        .toList();

                filteredAmberBAFs.addAll(regionAmberBAFs);
            }

            if(mConfig.CobaltDir != null)
            {
                List<CobaltRatio> regionCobaltRatios = mSampleData.CobaltRatios.stream()
                        .filter(x -> x.chromosome().equals(region.chromosome()) && region.start()<=x.position() && region.end()>=x.position())
                        .toList();

                filteredCobaltRatios.addAll(regionCobaltRatios);
            }

            if(mConfig.PurpleDir != null)
            {
                List<PurpleSegment> regionPurpleSegments = Lists.newArrayList();
                for(PurpleSegment purpleSegment : mSampleData.PurpleSegments)
                {
                    GenomeRegion purpleSegmentRegion = GenomeRegions.create(
                            purpleSegment.Chromosome, purpleSegment.PosStart, purpleSegment.PosEnd);

                    if(!region.overlaps(purpleSegmentRegion))
                        continue;

                    PurpleSegment newPurpleSegment = purpleSegment.withModifiedCoordinates(
                            purpleSegment.Chromosome,
                            Math.max(purpleSegment.PosStart, region.start()),
                            Math.min(purpleSegment.PosEnd, region.end())
                    );

                    regionPurpleSegments.add(newPurpleSegment);
                }

                filteredPurpleSegments.addAll(regionPurpleSegments);
            }
        }

        // Need to extend terminal segments past any current segments, links and exons and copy numbers
        List<VisSegment> segments = VisSegments.extendTerminals(
                0, filteredSegments, filteredLinks, positionsToCover, showSimpleSvSegments, mConfig.RefGenomeCoords);

        List<VisSvData> links = VisLinks.addFrame(segments, filteredLinks);

        ColorPicker color = colorPickerFactory.create(links);

        CircosData circosData = new CircosData(
                mCircosConfig,
                segments, links, copyNumbers, filteredExons, filteredFusions, filteredAmberBAFs, filteredCobaltRatios, filteredPurpleSegments,
                showSimpleSvSegments,  mConfig.IncludeFragileSites, mConfig.IncludeLineElements, fileId
        );

        CircosConfigWriter confWrite = new CircosConfigWriter(fileId, mConfig.OutputConfPath, circosData, mCircosConfig);
        FusionDataWriter fusionDataWriter = new FusionDataWriter(filteredFusions, filteredExons, filteredProteinDomains);

        mCallableConfigs.add(() -> new CircosDataWriter(color, fileId, mConfig.OutputConfPath, mCircosConfig, confWrite, circosData).write());
        if(!fusionDataWriter.finalExons().isEmpty())
        {
            mCallableConfigs.add(() -> fusionDataWriter.write(fileId, mConfig.OutputConfPath));
        }

        int minFrame = mCircosConfig.Step ? 0 : circosData.MaxFrame;
        for(int frame = minFrame; frame <= circosData.MaxFrame; frame++)
        {
            boolean plotFusion = !fusionDataWriter.finalExons().isEmpty();
            submitFrame(frame, plotFusion, circosData.labelSize(), fileId, confWrite);
        }
    }

    private void submitFrame(int frame, boolean hasFusion, double labelSize, String fileId, final CircosConfigWriter confWrite)
    {
        boolean plotFusion = hasFusion && !mConfig.Debug;
        boolean plotChromosome = !mConfig.Debug;

        mCallableConfigs.add(() -> confWrite.writeConfig(frame));
        mCallableImages.add(() -> createImageFrame(frame, labelSize, fileId, plotFusion, plotChromosome));
    }

    private Object createImageFrame(
            int frame, double labelSize, final String fileId, boolean plotFusion, boolean plotChromosome) throws Exception
    {
        String confFileName = fileId + ".circos." + String.format("%03d", frame) + ".conf";
        String outputFileName = fileId + "." + String.format("%03d", frame) + ".png";

        double rLabelSize = 1.2 * labelSize;

        Object circosResult = new CircosExecution(
                mConfig.CircosBin).generateCircos(mConfig.OutputConfPath + File.separator + confFileName,
                mConfig.OutputPlotPath, outputFileName);

        if(plotFusion)
        {
            int execution = new FusionExecution(fileId, outputFileName, mConfig.OutputConfPath, mConfig.OutputPlotPath)
                    .executeR(mCircosConfig, rLabelSize);

            if(execution != 0)
                throw new Exception("plotting error");
        }

        if(plotChromosome)
        {
            int execution = new ChromosomeRangeExecution(
                    fileId, mConfig.RefGenVersion, outputFileName, mConfig.OutputConfPath, mConfig.OutputPlotPath)
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

        try(SvVisualiser application = new SvVisualiser(configBuilder))
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
