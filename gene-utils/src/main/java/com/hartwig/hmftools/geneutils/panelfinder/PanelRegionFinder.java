package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.String.format;
import static java.lang.String.valueOf;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.getChromosomeFieldIndex;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_DEPTH_AVG;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_DEPTH_MAX;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_DEPTH_MIN;
import static com.hartwig.hmftools.common.region.HighDepthRegion.FLD_SAMPLE_COUNT;
import static com.hartwig.hmftools.common.region.TaggedRegion.loadRegionsFromBedFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
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
import static com.hartwig.hmftools.geneutils.panelfinder.PanelFinderConfig.DEFAULT_GENE_UPSTREAM_DISTANCE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.mappability.RegionQuality;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.HighDepthRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;
import com.hartwig.hmftools.common.sv.StartEndIterator;
import com.hartwig.hmftools.common.utils.StartEndPair;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;

import org.jetbrains.annotations.NotNull;

public class PanelRegionFinder
{
    private final PanelFinderConfig mConfig;

    private final Map<String,List<RegionData>> mChrRegions;
    private final Set<String> mPanelGeneNames;

    public PanelRegionFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new PanelFinderConfig(configBuilder);

        mChrRegions = Maps.newHashMap();
        mPanelGeneNames = Sets.newHashSet();
    }

    public void run()
    {
        GU_LOGGER.info("running panel region finder");

        // load high-depth regions
        loadHighDepthRegions();

        // load existing panel bed and merge
        loadPanelRegions();

        loadPanelGenes();

        // annotate with canonical transcript info
        annotateGeneExons();

        // annotate with genome mappability
        annotateMappability();

        applyFinalFilters();

        // write results
        writeResults();

        GU_LOGGER.info("panel region finder complete");
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
            int sampleCountIndex = fieldsIndexMap.get(FLD_SAMPLE_COUNT);
            int depthMinIndex = fieldsIndexMap.get(FLD_DEPTH_MIN);
            int depthMaxIndex = fieldsIndexMap.get(FLD_DEPTH_MAX);
            Integer depthAvgIndex = fieldsIndexMap.get(FLD_DEPTH_AVG);

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
                int depthAvg = depthAvgIndex != null ? Integer.parseInt(values[depthAvgIndex]) : 0;

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

                HighDepthRegion highDepthData = new HighDepthRegion(
                        chromosome, regionStart, regionEnd, depthMin, depthMax, depthAvg, sampleCount);

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

    private void loadPanelGenes()
    {
        if(mConfig.GeneIdFile != null)
        {
            mPanelGeneNames.addAll(loadDelimitedIdFile(mConfig.GeneIdFile, FLD_GENE_NAME, TSV_DELIM));
        }
        else if(mConfig.DriverGenePanel != null)
        {
            try
            {
                List<DriverGene> driverGenes = DriverGeneFile.read(mConfig.DriverGenePanel);
                driverGenes.forEach(x -> mPanelGeneNames.add(x.gene()));
            }
            catch(IOException e)
            {
                GU_LOGGER.error("invalid driver gene panel file: {}", e.toString());
                System.exit(1);
            }
        }

        if(!mPanelGeneNames.isEmpty())
        {
            GU_LOGGER.info("loaded {} panel gene names", mPanelGeneNames.size());
        }
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

            List<GeneTransData> geneTransDataList = Lists.newArrayListWithCapacity(geneDataList.size());

            // combine gene and transcript data
            for(GeneData geneData : geneDataList)
            {
                TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                if(transcriptData != null)
                {
                    geneTransDataList.add(new GeneTransData(geneData, transcriptData));
                }
            }

            List<GeneTransData> panelGeneDataList = geneTransDataList.stream()
                    .filter(x -> mPanelGeneNames.contains(x.geneData().GeneName)).collect(Collectors.toList());

            for(RegionData region : regions)
            {
                List<GeneTransData> overlappingTranscripts = geneTransDataList.stream()
                        .filter(x -> region.overlaps(x.geneData().Chromosome, x.transcriptData().TransStart, x.transcriptData().TransEnd))
                        .collect(Collectors.toList());

                GeneTransData panelGeneOverlap = panelGeneDataList.stream()
                        .filter(x -> regionOverlapsGeneBounds(region, x)).findFirst().orElse(null);

                if(panelGeneOverlap != null)
                    region.setPanelGene(panelGeneOverlap.geneData().GeneName);

                GeneData closestGene = null;
                ExonData closestExon = null;
                int closestNonOverlap = -1;

                if(!overlappingTranscripts.isEmpty())
                {
                    for(GeneTransData geneTransData : overlappingTranscripts)
                    {
                        GeneData geneData = geneTransData.geneData();
                        TranscriptData transcriptData = geneTransData.transcriptData();

                        // any region within exons will be annotated
                        for(int i = 0; i < transcriptData.exons().size(); ++i)
                        {
                            ExonData exon = transcriptData.exons().get(i);
                            ExonData nextExon = i < transcriptData.exons().size() - 1 ? transcriptData.exons().get(i + 1) : null;

                            if(positionsOverlap(region.start(), region.end(), exon.Start, exon.End))
                            {
                                region.addGeneExon(new GeneExonData(geneData.GeneName, exon.Rank, exon.Start, exon.End));
                            }
                            else if(region.geneExons().isEmpty() && nextExon != null
                            && region.start() > exon.End && region.end() < nextExon.Start)
                            {
                                int absDistance = min(abs(exon.End - region.start()), abs(nextExon.Start - region.end()));

                                if(closestNonOverlap == -1 || absDistance < closestNonOverlap)
                                {
                                    closestGene = geneData;
                                    closestExon = exon;
                                    closestNonOverlap = absDistance;
                                }
                            }
                        }
                    }
                }
                else
                {
                    for(GeneTransData geneTransData : geneTransDataList)
                    {
                        GeneData geneData = geneTransData.geneData();
                        TranscriptData transcriptData = geneTransData.transcriptData();

                        int absDistance = min(abs(transcriptData.TransEnd - region.start()), abs(transcriptData.TransStart - region.end()));

                        if(region.end() < geneData.GeneStart)
                        {
                            int permittedDistance = geneData.forwardStrand() ? mConfig.GeneUpstreamDistance : mConfig.GeneDownstreamDistance;

                            if(absDistance > permittedDistance)
                                break;
                        }
                        else
                        {
                            int permittedDistance = !geneData.forwardStrand() ? mConfig.GeneUpstreamDistance : mConfig.GeneDownstreamDistance;

                            if(absDistance > permittedDistance)
                                continue;
                        }

                        if(closestNonOverlap == -1 || absDistance < closestNonOverlap)
                        {
                            closestGene = geneData;
                            closestNonOverlap = absDistance;
                        }
                    }
                }

                // if no overlaps were found, then find the closest
                if(region.geneExons().isEmpty() && closestGene != null)
                {
                    String closeInfo;

                    if(closestExon != null)
                    {
                        closeInfo = format("%s exon(%d) distance(%d)", closestGene.GeneName, closestExon.Rank, closestNonOverlap);
                    }
                    else
                    {
                        boolean isUpstream = (region.end() < closestGene.GeneStart) == closestGene.forwardStrand();
                        closeInfo = format("%s %s distance(%d)", closestGene.GeneName, isUpstream ? "upstream" : "downstream", closestNonOverlap);
                    }

                    region.setClosestGeneInfo(closeInfo);
                }
            }
        }
    }

    private record GeneTransData(GeneData geneData, TranscriptData transcriptData) {}

    private boolean regionOverlapsGeneBounds(final RegionData region, final GeneTransData geneTransData)
    {
        int lowerGeneBounds = geneTransData.transcriptData().TransStart;
        int upperGeneBounds = geneTransData.transcriptData().TransEnd;

        if(geneTransData.geneData().forwardStrand())
        {
            lowerGeneBounds -= mConfig.GeneUpstreamDistance;
            upperGeneBounds += mConfig.GeneDownstreamDistance;
        }
        else
        {
            lowerGeneBounds -= mConfig.GeneDownstreamDistance;
            upperGeneBounds += mConfig.GeneUpstreamDistance;
        }

        return positionsOverlap(region.start(), region.end(), lowerGeneBounds, upperGeneBounds);
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

    private void applyFinalFilters()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<RegionData> regions = mChrRegions.get(chrStr);

            if(regions == null)
                continue;

            // should already be sorted but ensure
            Collections.sort(regions);

            int index = 0;

            while(index < regions.size())
            {
                RegionData region = regions.get(index);

                if(!region.panelRelated())
                {
                    if(mConfig.RequirePanelGene)
                    {
                        regions.remove(index);
                        continue;
                    }

                    // filter on mappability
                    double meanMappability = region.meanMappability();

                    if(mConfig.MinMappability > 0)
                    {
                        // apply the min mappability check if not in an existing panel region and not a known panel gene
                        if(meanMappability < mConfig.MinMappability)
                        {
                            regions.remove(index);
                            continue;
                        }
                    }
                }

                ++index;
            }

            if(mConfig.BackboneMinInterval > 0 || mConfig.BackboneMaxLength > 0)
            {
                index = 0;

                while(index < regions.size())
                {
                    RegionData region = regions.get(index);

                    if(!region.panelRelated())
                    {
                        boolean hasValidLength = mConfig.BackboneMaxLength == 0 || region.baseLength() <= mConfig.BackboneMaxLength;

                        boolean isCloseToOtherInterval = false;

                        if(mConfig.BackboneMinInterval > 0)
                        {
                            if(index > 0 && region.start() - regions.get(index - 1).end() < mConfig.BackboneMinInterval)
                            {
                                isCloseToOtherInterval = true;
                            }
                            else if(index < regions.size() - 1
                            && regions.get(index + 1).start() - region.end() < mConfig.BackboneMinInterval)
                            {
                                isCloseToOtherInterval = true;
                            }
                        }

                        if(!hasValidLength && isCloseToOtherInterval)
                        {
                            regions.remove(index);
                            continue;
                        }
                    }

                    ++index;
                }
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
            sj.add("HighDepthMax");
            sj.add("HighDepthAvg");
            sj.add("HighDepthSamples");
            sj.add("HighDepthInfo");

            sj.add("PanelRegionCount");
            sj.add("PanelRegionInfo");
            sj.add("PanelGene");
            sj.add("NonPanelRegionLength");
            sj.add("NearestPanelRegion");
            sj.add("NearestRegion");

            sj.add("GeneExonCount");
            sj.add("GeneExonInfo");
            sj.add("NearbyGeneInfo");

            sj.add("MappabilityAvg");
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

                for(int i = 0; i < regions.size(); ++i)
                {
                    RegionData region = regions.get(i);

                    sj = new StringJoiner(TSV_DELIM);
                    sj.add(region.Chromosome);
                    sj.add(valueOf(region.start()));
                    sj.add(valueOf(region.end()));

                    sj.add(valueOf(region.highDepths().size()));

                    int maxDepth = region.highDepths().stream().mapToInt(x -> x.DepthMax).max().orElse(0);
                    double avgDepth = region.highDepths().stream().mapToInt(x -> x.DepthAvg).average().orElse(0);
                    int maxSamples = region.highDepths().stream().mapToInt(x -> x.SampleCount).max().orElse(0);
                    sj.add(valueOf(maxDepth));
                    sj.add(format("%.0f", avgDepth));
                    sj.add(valueOf(maxSamples));
                    sj.add(valueOf(RegionData.toString(region.highDepths())));

                    sj.add(valueOf(region.panelRegions().size()));
                    sj.add(valueOf(PanelData.toString(region.panelRegions())));

                    // report closest panel region if not one
                    int closestPanelRegion = -1;
                    int closestRegion = -1;

                    // search up and down from the current region
                    for(int j = 0; j <= 1; ++j)
                    {
                        boolean searchDown = (j == 0);
                        int nextIndex = searchDown ? i - 1 : i + 1;

                        while(nextIndex >= 0 && nextIndex < regions.size())
                        {
                            RegionData nextRegion = regions.get(nextIndex);
                            int distance = searchDown ? region.start() - nextRegion.end() : nextRegion.start() - region.end();

                            if(closestRegion < 0 || distance < closestRegion)
                                closestRegion = distance;

                            if(!nextRegion.panelRegions().isEmpty())
                            {
                                 if(closestPanelRegion < 0 || distance < closestPanelRegion)
                                     closestPanelRegion = distance;

                                break;
                            }

                            nextIndex += searchDown ?  -1 : 1;
                        }
                    }

                    sj.add(region.panelGeneName());

                    // length of the region not covered by panel regions
                    int regionLength = region.baseLength();
                    int panelRegionsLength = region.panelRegions().stream().mapToInt(x -> x.baseLength()).sum();
                    int nonPanelLength = max(regionLength - panelRegionsLength, 0);
                    sj.add(valueOf(nonPanelLength));

                    sj.add(valueOf(closestPanelRegion));
                    sj.add(valueOf(closestRegion));

                    sj.add(valueOf(region.geneExons().size()));
                    sj.add(valueOf(GeneExonData.toString(region.geneExons())));

                    sj.add(valueOf(region.closestGeneInfo()));

                    double minMappability = region.mappabilityScores().stream().mapToDouble(x -> x.Quality).min().orElse(0);
                    double maxMappability = region.mappabilityScores().stream().mapToDouble(x -> x.Quality).max().orElse(0);
                    double meanMappability = region.meanMappability();

                    sj.add(format("%.3f", meanMappability));
                    sj.add(format("%.3f", minMappability));
                    sj.add(format("%.3f", maxMappability));

                    writer.write(sj.toString());
                    writer.newLine();

                    if(bedWriter != null)
                    {
                        sj = new StringJoiner(TSV_DELIM);
                        sj.add(region.Chromosome);
                        sj.add(valueOf(region.start() - 1)); // since as BED
                        sj.add(valueOf(region.end()));
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
