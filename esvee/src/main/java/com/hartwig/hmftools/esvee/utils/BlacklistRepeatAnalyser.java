package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig.addGenePanelOption;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.KNOWN_FUSIONS_FILE;
import static com.hartwig.hmftools.common.fusion.KnownFusionCache.addKnownFusionFileOption;
import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.sv.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.immune.ImmuneRegions.getIgRegion;
import static com.hartwig.hmftools.common.mappability.UnmappedRegions.UNMAP_REGIONS_FILE;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.PARTITION_SIZE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_CHR_PARTITION_SIZE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.sv.RepeatMaskData;
import com.hartwig.hmftools.common.mappability.UnmappedRegions;
import com.hartwig.hmftools.common.mappability.UnmappingRegion;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.prep.BlacklistLocations;

import org.jetbrains.annotations.NotNull;

public class BlacklistRepeatAnalyser
{
    private final String mOutputFile;

    private final BlacklistLocations mBlacklistLocations;
    private final RepeatMaskAnnotations mRepeatMaskAnnotations;
    private final Map<String, List<UnmappingRegion>> mUnmapRegionsMap;
    private final BufferedWriter mWriter;
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenVersion;

    private final int mBaseWindowLength;
    private final double mMinDominantPercent;

    private final boolean mRemoveGeneOverlaps;
    private final boolean mWriteBlacklistBed;
    private final Map<String, List<GeneInfo>> mChrGenicRegions;
    private final List<ChrBaseRegion> mFinalBlacklistRegions;

    private final int mPartitionSize;
    private final int mThreads;
    private final SpecificRegions mSpecificChrRegions;

    public static final String BLACKLIST_BED = "blacklist_bed";

    private static final String OUTPUT_FILE = "output_file";
    private static final String MIN_DOMINANT_PERC = "min_dominant_perc";
    protected static final String BASE_WINDOW_LENGTH = "base_window_length";
    private static final String REMOVE_GENE_OVERLAPS = "remove_gene_overlaps";
    private static final String WRITE_BLACKLIST_BED = "write_blacklist_bed";

    private static final int DEFAULT_BASE_WINDOW_LENGTH = 100;
    private static final double DEFAULT_MIN_DOMINANT_PERC = 0.95;

    public BlacklistRepeatAnalyser(final ConfigBuilder configBuilder)
    {
        mBlacklistLocations = new BlacklistLocations(configBuilder.getValue(BLACKLIST_BED));
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mWriter = initialiseWriter();

        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
        mRefGenVersion = deriveRefGenomeVersion(mRefGenome);

        mRepeatMaskAnnotations = new RepeatMaskAnnotations();
        if(configBuilder.hasValue(REPEAT_MASK_FILE))
        {
            if(!mRepeatMaskAnnotations.load(configBuilder.getValue(REPEAT_MASK_FILE), mRefGenVersion))
                System.exit(1);
        }

        if(configBuilder.hasValue(UNMAP_REGIONS_FILE))
        {
            mUnmapRegionsMap = UnmappedRegions.loadUnmapRegions(configBuilder.getValue(UNMAP_REGIONS_FILE));
        }
        else
        {
            mUnmapRegionsMap = Collections.emptyMap();
        }

        mChrGenicRegions = Maps.newHashMap();
        loadGenicInfo(configBuilder);

        mFinalBlacklistRegions = Lists.newArrayList();

        mPartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        mBaseWindowLength = configBuilder.getInteger(BASE_WINDOW_LENGTH);
        mMinDominantPercent = configBuilder.getDecimal(MIN_DOMINANT_PERC);
        mRemoveGeneOverlaps = configBuilder.hasFlag(REMOVE_GENE_OVERLAPS);
        mWriteBlacklistBed = configBuilder.hasFlag(WRITE_BLACKLIST_BED);

        mThreads = parseThreads(configBuilder);
        mSpecificChrRegions = SpecificRegions.from(configBuilder);
    }

    public void run()
    {
        SV_LOGGER.info("analysing blacklist regions");

        long startTimeMs = System.currentTimeMillis();

        RefGenomeCoordinates coordinates = mRefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        List<ChrBaseRegion> regions = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenVersion.versionedChromosome(chromosome.toString());

            if(mSpecificChrRegions.excludeChromosome(chrStr))
                continue;

            List<ChrBaseRegion> chrRegions = buildPartitions(chrStr, coordinates.length(chrStr), mPartitionSize);

            for(ChrBaseRegion region : chrRegions)
            {
                if(mSpecificChrRegions.hasFilters() && !mSpecificChrRegions.includeRegion(region))
                    continue;

                regions.add(region);
            }
        }

        List<RegionTask> regionTasks = regions.stream().map(x -> new RegionTask(x)).collect(Collectors.toList());

        List<Callable<Void>> callableList = regionTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        closeBufferedWriter(mWriter);

        if(mWriteBlacklistBed)
        {
            ChrBaseRegion.checkMergeOverlaps(mFinalBlacklistRegions, true);

            try
            {
                String outputBed = mOutputFile.replaceAll(TSV_EXTENSION, ".bed");

                BufferedWriter writer = createBufferedWriter(outputBed, false);

                StringJoiner sj = new StringJoiner(TSV_DELIM);

                sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END);

                writer.write(sj.toString());
                writer.newLine();

                for(ChrBaseRegion region : mFinalBlacklistRegions)
                {
                    StringJoiner regionData = new StringJoiner(TSV_DELIM);
                    regionData.add(region.Chromosome);
                    regionData.add(String.valueOf(region.start() - 1));  // write as a BED file, so note the -1 on the start
                    regionData.add(String.valueOf(region.end()));
                    writer.write(regionData.toString());
                    writer.newLine();
                }

                writer.close();
            }
            catch(IOException e)
            {
                SV_LOGGER.error(" failed to write final BED file: {}", e.toString());
                System.exit(1);
            }
        }

        SV_LOGGER.info("Blacklist region analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private class RegionTask implements Callable<Void>
    {
        private final ChrBaseRegion mRegion;

        private final List<RepeatMaskData> mRepeatMaskRegions;
        private final List<BaseRegion> mBlackListRegions;
        private final List<UnmappingRegion> mUnmapRegions;
        private final List<GeneInfo> mGeneRegions;

        public RegionTask(final ChrBaseRegion region)
        {
            mRegion = region;
            mRepeatMaskRegions = mRepeatMaskAnnotations.findMatches(region);

            mBlackListRegions = Lists.newArrayList();
            List<BaseRegion> chrBlacklistRegions = mBlacklistLocations.getRegions(region.Chromosome);

            if(chrBlacklistRegions != null)
            {
                chrBlacklistRegions.stream()
                        .filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                        .forEach(x -> mBlackListRegions.add(x));
            }

            mUnmapRegions = Lists.newArrayList();
            List<UnmappingRegion> chrUnmapRegions = mUnmapRegionsMap.get(region.Chromosome);

            if(chrUnmapRegions != null)
            {
                chrUnmapRegions.stream()
                        .filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                        .forEach(x -> mUnmapRegions.add(x));
            }

            mGeneRegions = Lists.newArrayList();
            List<GeneInfo> chrGeneRegions = mChrGenicRegions.get(region.Chromosome);

            if(chrGeneRegions != null)
            {
                chrGeneRegions.stream()
                        .filter(x -> positionsOverlap(mRegion.start(), mRegion.end(), x.start(), x.end()))
                        .forEach(x -> mGeneRegions.add(x));
            }
        }

        @Override
        public Void call()
        {
            int currentWindowStart = mRegion.start();
            int windowEnd = currentWindowStart + mBaseWindowLength - 1;

            int minRegionEnd = mRegion.end() - (int) (mBaseWindowLength * 0.1);

            while(windowEnd <= minRegionEnd)
            {
                analyseWindow(new BaseRegion(currentWindowStart, min(windowEnd, mRegion.end())));

                currentWindowStart = windowEnd + 1;
                windowEnd = currentWindowStart + mBaseWindowLength - 1;
            }

            return null;
        }

        private void analyseWindow(final BaseRegion region)
        {
            // ignore regions already covered by an unmapping region
            for(UnmappingRegion unmapRegion : mUnmapRegions)
            {
                // check for same required overlap as used in Redux to skip a read
                if(hasSufficientOverlap(region, unmapRegion, 0.9))
                    continue;
            }

            byte[] refBases = mRefGenome.getBases(mRegion.Chromosome, region.start(), region.end());

            short[] baseCounts = new short[Nucleotides.DNA_BASES.length];

            for(int i = 0; i < refBases.length; ++i)
            {
                int baseIndex = Nucleotides.baseIndex(refBases[i]);

                if(baseIndex >= 0 && baseIndex < baseCounts.length)
                    ++baseCounts[baseIndex];
            }

            char topBase = 0;
            char secondBase = 0;
            float topBaseCount = 0;
            float secondBaseCount = 0;

            for(int i = 0; i < baseCounts.length; ++i)
            {
                if(baseCounts[i] > topBaseCount)
                {
                    if(topBaseCount > secondBaseCount)
                    {
                        secondBaseCount = topBaseCount;
                        secondBase = topBase;
                    }

                    topBaseCount = baseCounts[i];
                    topBase = Nucleotides.DNA_BASES[i];
                }
                else if(baseCounts[i] == topBaseCount || baseCounts[i] > secondBaseCount)
                {
                    secondBaseCount = baseCounts[i];
                    secondBase = Nucleotides.DNA_BASES[i];
                }
            }

            String dominantBases;
            double dominantBasePercent;

            double topBasePerc = topBaseCount / (double) mBaseWindowLength;
            double topAndSecondBasePerc = (topBaseCount + secondBaseCount) / (double) mBaseWindowLength;

            if(topBasePerc >= mMinDominantPercent)
            {
                dominantBases = String.valueOf(topBase);
                dominantBasePercent = topBasePerc;
            }
            else if(topAndSecondBasePerc >= mMinDominantPercent)
            {
                dominantBases = format("%c_%c", topBase, secondBase);
                dominantBasePercent = topAndSecondBasePerc;
            }
            else
            {
                return;
            }

            List<RepeatMaskData> repeatMasks = mRepeatMaskRegions.stream()
                    .filter(x -> hasSufficientOverlap(region, x.Region, 0.5)).collect(Collectors.toList());

            List<BaseRegion> blacklistRegions = mBlackListRegions.stream()
                    .filter(x -> hasSufficientOverlap(region, x, 0.5)).collect(Collectors.toList());

            List<GeneInfo> genicOverlaps = mGeneRegions.stream()
                    .filter(x -> positionsOverlap(x.start(), x.end(), region.start(), region.end())).collect(Collectors.toList());

            writeRegionData(mRegion.Chromosome, region, dominantBases, dominantBasePercent, repeatMasks, blacklistRegions, genicOverlaps);

            if(mWriteBlacklistBed)
            {
                if(!mRemoveGeneOverlaps || genicOverlaps.isEmpty())
                {
                    addBlacklistRegion(new ChrBaseRegion(mRegion.Chromosome, region.start(), region.end()));
                }
            }
        }
    }

    private synchronized void addBlacklistRegion(final ChrBaseRegion region)
    {
        mFinalBlacklistRegions.add(region);
    }

    private boolean hasSufficientOverlap(final BaseRegion window, final BaseRegion otherRegion, final double requiredPerc)
    {
        if(!otherRegion.overlaps(window))
            return false;

        int overlap = min(otherRegion.end(), window.end()) - max(otherRegion.start(), window.start()) + 1;
        return overlap >= requiredPerc * mBaseWindowLength;
    }

    private void loadGenicInfo(final ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(ENSEMBL_DATA_DIR))
            return;

        List<DriverGene> driverGenes = DriverGenePanelConfig.loadDriverGenes(configBuilder);
        KnownFusionCache knownFusionCache = new KnownFusionCache();
        knownFusionCache.loadFile(configBuilder.getValue(KNOWN_FUSIONS_FILE));

        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(configBuilder.getValue(ENSEMBL_DATA_DIR), mRefGenVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        Set<String> addedGenes = driverGenes.stream().map(x -> x.gene()).collect(Collectors.toSet());

        List<KnownFusionData> knownFusionData = knownFusionCache.getDataByType(KNOWN_PAIR);

        if(knownFusionData != null)
        {
            knownFusionData.forEach(x -> addedGenes.add(x.FiveGene));
            knownFusionData.forEach(x -> addedGenes.add(x.ThreeGene));
        }

        for(String geneName : addedGenes)
        {
            GeneData geneData = ensemblDataCache.getGeneDataByName(geneName);

            if(geneData == null)
                continue;

            TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

            if(transcriptData == null)
                continue;

            List<GeneInfo> geneRegions = mChrGenicRegions.get(geneData.Chromosome);

            if(geneRegions == null)
            {
                geneRegions = Lists.newArrayList();
                mChrGenicRegions.put(geneData.Chromosome, geneRegions);
            }

            if(transcriptData.posStrand())
            {
                geneRegions.add(new GeneInfo(transcriptData.
                        TransStart - 10000, transcriptData.TransStart - 1, format("%s_UPSTREAM", geneName)));
            }
            else
            {
                geneRegions.add(new GeneInfo(transcriptData.
                        TransEnd + 1, transcriptData.TransEnd + 10000, format("%s_UPSTREAM", geneName)));
            }

            for(int i = 0; i < transcriptData.exons().size(); ++i)
            {
                ExonData exonData = transcriptData.exons().get(i);

                geneRegions.add(new GeneInfo(exonData.Start, exonData.End, format("%s_EXON_%d", geneName, exonData.Rank)));

                if(i < transcriptData.exons().size() - 1)
                {
                    ExonData nextExon = transcriptData.exons().get(i + 1);

                    geneRegions.add(new GeneInfo(
                            exonData.End + 1, nextExon.Start - 1, format("%s_INTRON_%d", geneName, exonData.Rank)));
                }
            }
        }

        // manually create genic regions for IG
        List<String> igGenes = List.of("IGH", "IGK", "IGL");

        for(String geneName : igGenes)
        {
            ChrBaseRegion igRegion = getIgRegion(geneName, mRefGenVersion);

            List<GeneInfo> geneRegions = mChrGenicRegions.get(igRegion.Chromosome);

            if(geneRegions == null)
            {
                geneRegions = Lists.newArrayList();
                mChrGenicRegions.put(igRegion.Chromosome, geneRegions);
            }

            geneRegions.add(new GeneInfo(igRegion.start(), igRegion.end(), geneName));
        }
    }

    private class GeneInfo extends BaseRegion
    {
        public final String Info;

        public GeneInfo(final int posStart, final int posEnd, final String info)
        {
            super(posStart, posEnd);
            Info = info;
        }

        @Override public String toString() { return format("%s:%d-%d", Info, start(), end()); }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END);

            sj.add("DominantBases").add("Concentration").add("BlacklistRegions").add("RepeatInfo").add("GeneInfo");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeRegionData(
            final String chromosome, final BaseRegion region, final String dominantBases, final double basePercent,
            final List<RepeatMaskData> repeatMasks, final List<BaseRegion> blacklistRegions, final List<GeneInfo> genicOverlaps)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(chromosome).add(String.valueOf(region.start())).add(String.valueOf(region.end()));

            sj.add(dominantBases).add(format("%.2f", basePercent));

            String blacklistInfo = blacklistRegions.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));

            String repeatInfo = repeatMasks.stream()
                    .map(x -> format("%s:%d-%d", x.ClassType, x.Region.start(), x.Region.end())).collect(Collectors.joining(ITEM_DELIM));

            String geneOverlaps = genicOverlaps.stream()
                    .map(x -> format("%s_%d:%d", x.Info, x.start(), x.end())).collect(Collectors.joining(ITEM_DELIM));

            sj.add(blacklistInfo).add(repeatInfo).add(geneOverlaps);

            mWriter.write(sj.toString());
            mWriter.newLine();

        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addRefGenomeFile(configBuilder, true);
        configBuilder.addPath(BLACKLIST_BED, false, "Blacklist BED file");
        UnmappedRegions.registerConfig(configBuilder);
        RepeatMaskAnnotations.addConfig(configBuilder);

        addGenePanelOption(configBuilder, false);
        addKnownFusionFileOption(configBuilder);
        addEnsemblDir(configBuilder);

        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addInteger(PARTITION_SIZE, "Read length", DEFAULT_CHR_PARTITION_SIZE);

        configBuilder.addInteger(BASE_WINDOW_LENGTH, "Base window length for analysis", DEFAULT_BASE_WINDOW_LENGTH);
        configBuilder.addDecimal(MIN_DOMINANT_PERC, "Min dominant concerntration of 1 or 2 bases", DEFAULT_MIN_DOMINANT_PERC);
        configBuilder.addFlag(REMOVE_GENE_OVERLAPS, "Remove regions which overlap a genic region");
        configBuilder.addFlag(WRITE_BLACKLIST_BED, "Write blacklist BED file for use in Esvee Prep");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BlacklistRepeatAnalyser blacklistRepeatAnalyser = new BlacklistRepeatAnalyser(configBuilder);
        blacklistRepeatAnalyser.run();
    }
}
