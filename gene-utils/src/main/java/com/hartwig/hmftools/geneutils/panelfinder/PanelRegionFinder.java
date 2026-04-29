package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.TaggedRegion.loadRegionsFromBedFile;
import static com.hartwig.hmftools.common.utils.Doubles.median;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
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

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
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
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.mappability.RegionQuality;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class PanelRegionFinder
{
    private static final String HIGH_DEPTH_FILE = "high_depth_file";
    private static final String OUTPUT_FILE = "output_file";

    private final String mHighDepthFile;
    private final String mTargetRegionsBed;
    private final String mEnsemblDataPath;
    private final String mMappabilityProfileFile;
    private final RefGenomeVersion mRefGenomeVersion;
    private final String mOutputFile;

    private final Map<String,List<RegionData>> mChrRegions;

    public PanelRegionFinder(final ConfigBuilder configBuilder)
    {
        mHighDepthFile = configBuilder.getValue(HIGH_DEPTH_FILE);
        mTargetRegionsBed = configBuilder.getValue(TARGET_REGIONS_BED);
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mRefGenomeVersion = RefGenomeVersion.from(configBuilder);

        mEnsemblDataPath = configBuilder.getValue(ENSEMBL_DATA_DIR);
        mMappabilityProfileFile = configBuilder.getValue(ProbeQualityProfile.CFG_PROBE_QUALITY_FILE);

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
            List<String> lines = Files.readAllLines(Paths.get(mHighDepthFile));

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

            int count = 0;

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                HighDepthData highDepthData = new HighDepthData(
                        new ChrBaseRegion(values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex])),
                        Integer.parseInt(values[sampleCountIndex]),
                        Integer.parseInt(values[depthMinIndex]), Integer.parseInt(values[depthMaxIndex]));

                ++count;

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

            GU_LOGGER.info("loaded {} high-depth regions from file({})", count, mHighDepthFile);
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to load high-depth regions file({}): {}", mHighDepthFile, e.toString());
            System.exit(1);
        }
    }

    private void loadPanelRegions()
    {
        if(mTargetRegionsBed == null)
            return;

        Map<Chromosome,List<TaggedRegion>> chrPanelRegions = loadRegionsFromBedFile(mTargetRegionsBed);

        GU_LOGGER.info("loaded {} panel regions from file({})",
                chrPanelRegions.values().stream().mapToInt(x -> x.size()).sum(), mTargetRegionsBed);

        for(Map.Entry<Chromosome,List<TaggedRegion>> entry : chrPanelRegions.entrySet())
        {
            String chromosome = mRefGenomeVersion.versionedChromosome(entry.getKey().toString());

            List<RegionData> regions = mChrRegions.get(chromosome);
            Collections.sort(regions);

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
        if(mEnsemblDataPath == null)
            return;

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(mEnsemblDataPath, mRefGenomeVersion);

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
        if(mMappabilityProfileFile == null)
            return;

        ProbeQualityProfile probeQualityProfile = ProbeQualityProfile.loadFromResourceFile(mMappabilityProfileFile);

        GU_LOGGER.debug("loaded genome-mappability file({})", mMappabilityProfileFile);

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
        GU_LOGGER.info("writing {} regions to file({})", mChrRegions.values().stream().mapToInt(x -> x.size()).sum(), mOutputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile);

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
                String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

                List<RegionData> regions = mChrRegions.get(chrStr);

                if(regions == null)
                    continue;

                // should already be sorted but ensure
                Collections.sort(regions);

                for(RegionData region : regions)
                {
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
                    sj.add(format("%.3f", region.meanMappability()));
                    sj.add(format("%.3f", minMappability));
                    sj.add(format("%.3f", maxMappability));

                    writer.write(sj.toString());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write output file({}): {}", mOutputFile, e.toString());
            System.exit(1);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(HIGH_DEPTH_FILE, true, "Input regions TSV file - header is optional");
        configBuilder.addPath(TARGET_REGIONS_BED, false, TARGET_REGIONS_BED_DESC);
        EnsemblDataCache.addEnsemblDir(configBuilder, false);
        addRefGenomeVersion(configBuilder);
        configBuilder.addConfigItem(OUTPUT_FILE, false, "Output filename");
        ProbeQualityProfile.registerConfig(configBuilder);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PanelRegionFinder panelRegionFinder = new PanelRegionFinder(configBuilder);
        panelRegionFinder.run();
    }
}
