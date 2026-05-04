package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.TaggedRegion.loadRegionsFromBedFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.panelfinder.PanelFinderConfig.CHROMOSOME_Y_SAMPLE_FRACTION;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.mappability.RegionQuality;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class PanelRegionFinder
{
    private final PanelFinderConfig mConfig;

    private final Map<String,List<RegionData>> mChrRegions;

    public PanelRegionFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelFinderConfig(configBuilder);

        mChrRegions = Maps.newHashMap();
    }

    public void run()
    {
        GU_LOGGER.info("running panel region finder");

        // load high-depth regions
        loadHighDepthRegions();

        // load existing panel bed and merge
        loadPanelRegions();

        // annotate with canonical transcript info
        annotateGeneExons();

        // annotate with genome mappability
        annotateMappability();

        // write results
        writeResults();

        GU_LOGGER.info("panel region finder complete");
    }

    private static boolean isBedFile(final String filename)
    {
        return filename.endsWith(".bed") || filename.endsWith(".bed.gz");
    }

    private void loadHighDepthRegions()
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(mConfig.HighDepthFile));

            // check for headers
            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_POS_END);
            int sampleCountIndex = fieldsIndexMap.get("SampleCount");
            int depthMinIndex = fieldsIndexMap.get("DepthMin");
            int depthMaxIndex = fieldsIndexMap.get("DepthMax");

            double minSamplesChromosomeY = CHROMOSOME_Y_SAMPLE_FRACTION * mConfig.MinSampleCount;

            int count = 0;
            int filtered = 0;

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String chromosome = values[chrIndex];
                int regionStart = Integer.parseInt(values[posStartIndex]);
                int regionEnd = Integer.parseInt(values[posEndIndex]);
                int sampleCount = Integer.parseInt(values[sampleCountIndex]);
                int depthMin = Integer.parseInt(values[depthMinIndex]);
                int depthMax = Integer.parseInt(values[depthMaxIndex]);

                ++count;

                if(mConfig.MinSampleCount > 0)
                {
                    if(HumanChromosome.fromString(chromosome) == _Y)
                    {
                        if(sampleCount < minSamplesChromosomeY)
                        {
                            ++filtered;
                            continue;
                        }
                    }
                    else if(sampleCount < mConfig.MinSampleCount)
                    {
                        ++filtered;
                        continue;
                    }
                }

                if(mConfig.HighDepthTrimCount > 0)
                {
                    regionStart += mConfig.HighDepthTrimCount;
                    regionEnd -= mConfig.HighDepthTrimCount;

                    if(regionStart > regionEnd)
                    {
                        ++filtered;
                        continue;
                    }
                }

                HighDepthData highDepthData = new HighDepthData(
                        new ChrBaseRegion(chromosome, regionStart, regionEnd), sampleCount, depthMin, depthMax);

                List<RegionData> regions = mChrRegions.get(highDepthData.Chromosome);

                if(regions == null)
                {
                    regions = Lists.newArrayList();
                    mChrRegions.put(highDepthData.Chromosome, regions);
                }

                RegionData regionData = new RegionData(highDepthData);
                regionData.addHighDepth(highDepthData);
                regions.add(regionData);
            }

            GU_LOGGER.info("loaded {} high-depth regions, filtered({}), from file({})", count, filtered, mConfig.HighDepthFile);
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to load high-depth regions file({}): {}", mConfig.HighDepthFile, e.toString());
            System.exit(1);
        }
    }

    private void loadPanelRegions()
    {
        if(mConfig.TargetRegionsBed == null)
            return;

        Map<Chromosome,List<TaggedRegion>> chrPanelRegions = loadRegionsFromBedFile(mConfig.TargetRegionsBed);

        GU_LOGGER.info("loaded {} panel regions from file({})",
                chrPanelRegions.values().stream().mapToInt(x -> x.size()).sum(), mConfig.TargetRegionsBed);

        for(Map.Entry<Chromosome,List<TaggedRegion>> entry : chrPanelRegions.entrySet())
        {
            String chromosome = mConfig.RefGenVersion.versionedChromosome(entry.getKey().toString());
            List<TaggedRegion> panelRegions = entry.getValue();
            Collections.sort(panelRegions);

            List<RegionData> regions = mChrRegions.get(chromosome);

            if(regions == null)
            {
                regions = Lists.newArrayList();
                mChrRegions.put(chromosome, regions);
            }

            for(TaggedRegion taggedRegion : entry.getValue())
            {
                PanelData panelData = new PanelData(taggedRegion, taggedRegion.Tag);
                mergePanelRegion(panelData, regions);
            }
        }
    }

    private void mergePanelRegion(final PanelData panelData, final List<RegionData> regions)
    {
        int index = 0;

        while(index < regions.size())
        {
            RegionData region = regions.get(index);

            if(panelData.start() > region.end())
            {
                ++index;
                continue;
            }

            if(panelData.end() < region.start())
            {
                RegionData regionData = new RegionData(panelData);
                regionData.addPanelData(panelData);
                regions.add(index, regionData);
                return;
            }

            // otherwise merge
            region.addPanelData(panelData);

            // remove any subsequent regions now covered
            int nextIndex = index + 1;
            while(nextIndex < regions.size())
            {
                RegionData nextRegion = regions.get(nextIndex);

                if(nextRegion.start() <= region.end())
                {
                    region.mergeRegion(nextRegion);
                    regions.remove(nextIndex);
                }
                else
                {
                    break;
                }
            }

            return;
        }

        RegionData regionData = new RegionData(panelData);
        regionData.addPanelData(panelData);
        regions.add(regionData);
    }

    private void annotateGeneExons()
    {
        if(mConfig.EnsemblDataPath == null)
            return;

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(mConfig.EnsemblDataPath, mConfig.RefGenVersion);

        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        for(Map.Entry<String,List<RegionData>> entry : mChrRegions.entrySet())
        {
            List<RegionData> regions = entry.getValue();

            List<GeneData> geneDataList = ensemblDataCache.getChrGeneDataMap().get(entry.getKey());

            if(geneDataList == null)
                continue;

            for(RegionData region : regions)
            {
                List<GeneData> overlappingGenes = geneDataList.stream()
                        .filter(x -> region.overlaps(x.Chromosome, x.GeneStart, x.GeneEnd)).collect(Collectors.toList());

                if(overlappingGenes.isEmpty())
                    continue;

                for(GeneData geneData : overlappingGenes)
                {
                    TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                    if(transcriptData == null)
                        continue;

                    for(ExonData exon : transcriptData.exons())
                    {
                        if(positionsOverlap(region.start(), region.end(), exon.Start, exon.End))
                            region.addGeneExon(new GeneExonData(geneData.GeneName, exon.Rank, exon.Start, exon.End));
                    }
                }
            }
        }
    }

    private void annotateMappability()
    {
        if(mConfig.MappabilityProfileFile == null)
            return;

        ProbeQualityProfile probeQualityProfile = ProbeQualityProfile.loadFromResourceFile(mConfig.MappabilityProfileFile);

        GU_LOGGER.debug("loaded genome-mappability file({})", mConfig.MappabilityProfileFile);

        for(Map.Entry<String,List<RegionData>> entry : mChrRegions.entrySet())
        {
            List<RegionData> regions = entry.getValue();

            for(RegionData region : regions)
            {
                List<RegionQuality> regionQualities = probeQualityProfile.findRegionQualities(region);
                region.mappabilityScores().addAll(regionQualities);
            }
        }
    }

    private void writeResults()
    {
        GU_LOGGER.info("writing {} regions to file({})", mChrRegions.values().stream().mapToInt(x -> x.size()).sum(), mConfig.OutputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.OutputFile);
            BufferedWriter bedWriter = mConfig.OutputBed != null ? createBufferedWriter(mConfig.OutputBed) : null;

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME);
            sj.add(FLD_REGION_START);
            sj.add(FLD_REGION_END);

            sj.add("HighDepthCount");
            sj.add("HighDepthMaxDepth");
            sj.add("HighDepthMaxSamples");
            sj.add("HighDepthInfo");

            sj.add("PanelRegionCount");
            sj.add("PanelRegionInfo");

            sj.add("GeneExonCount");
            sj.add("GeneExonInfo");

            sj.add("MappabilityMedian");
            sj.add("MappabilityMin");
            sj.add("MappabilityMax");

            writer.write(sj.toString());
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

                List<RegionData> regions = mChrRegions.get(chrStr);

                if(regions == null)
                    continue;

                // should already be sorted but ensure
                Collections.sort(regions);

                for(RegionData region : regions)
                {
                    // filter on mappability
                    double meanMappability = region.meanMappability();

                    if(mConfig.MinMappability > 0)
                    {
                        if(region.panelRegions().isEmpty() && meanMappability < mConfig.MinMappability)
                            continue;
                    }

                    sj = new StringJoiner(TSV_DELIM);
                    sj.add(region.Chromosome);
                    sj.add(String.valueOf(region.start()));
                    sj.add(String.valueOf(region.end()));

                    sj.add(String.valueOf(region.highDepths().size()));

                    int maxDepth = region.highDepths().stream().mapToInt(x -> x.MaxDepth).max().orElse(0);
                    int maxSamples = region.highDepths().stream().mapToInt(x -> x.SampleCount).max().orElse(0);
                    sj.add(String.valueOf(maxDepth));
                    sj.add(String.valueOf(maxSamples));
                    sj.add(String.valueOf(HighDepthData.toString(region.highDepths())));

                    sj.add(String.valueOf(region.panelRegions().size()));
                    sj.add(String.valueOf(PanelData.toString(region.panelRegions())));

                    sj.add(String.valueOf(region.geneExons().size()));
                    sj.add(String.valueOf(GeneExonData.toString(region.geneExons())));

                    double minMappability = region.mappabilityScores().stream().mapToDouble(x -> x.Quality).min().orElse(0);
                    double maxMappability = region.mappabilityScores().stream().mapToDouble(x -> x.Quality).max().orElse(0);
                    sj.add(format("%.3f", meanMappability));
                    sj.add(format("%.3f", minMappability));
                    sj.add(format("%.3f", maxMappability));

                    writer.write(sj.toString());
                    writer.newLine();

                    if(bedWriter != null)
                    {
                        sj = new StringJoiner(TSV_DELIM);
                        sj.add(region.Chromosome);
                        sj.add(String.valueOf(region.start() - 1)); // since as BED
                        sj.add(String.valueOf(region.end()));
                        sj.add(region.label());

                        bedWriter.write(sj.toString());
                        bedWriter.newLine();
                    }
                }
            }

            writer.close();

            if(bedWriter != null)
                bedWriter.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write output file({}): {}", mConfig.OutputFile, e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        PanelFinderConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PanelRegionFinder panelRegionFinder = new PanelRegionFinder(configBuilder);
        panelRegionFinder.run();
    }
}
