package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations.REPEAT_MASK_FILE;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.region.UnmappedRegions.UNMAP_REGIONS_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.BLACKLIST_BED;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.PARTITION_SIZE;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_CHR_PARTITION_SIZE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.gripss.RepeatMaskAnnotations;
import com.hartwig.hmftools.common.gripss.RepeatMaskData;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.region.UnmappedRegions;
import com.hartwig.hmftools.common.region.UnmappingRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.prep.BlacklistLocations;

import org.jetbrains.annotations.NotNull;

public class BlacklistRepeatAnalyser
{
    private final String mOutputFile;

    private final BlacklistLocations mBlacklistLocations;
    private final RepeatMaskAnnotations mRepeatMaskAnnotations;
    private final Map<String,List<UnmappingRegion>> mUnmapRegionsMap;
    private final BufferedWriter mWriter;
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenVersion;

    private final int mBaseWindowLength;
    private final double mMinDominantPercent;

    private final int mPartitionSize;
    private final int mThreads;
    private final SpecificRegions mSpecificChrRegions;

    private static final String OUTPUT_FILE = "output_file";
    private static final String MIN_DOMINANT_PERC = "min_dominant_perc";
    private static final String BASE_WINDOW_LENGTH = "base_window_length";

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

        mPartitionSize = configBuilder.getInteger(PARTITION_SIZE);
        mBaseWindowLength = configBuilder.getInteger(BASE_WINDOW_LENGTH);
        mMinDominantPercent = configBuilder.getDecimal(MIN_DOMINANT_PERC);

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

        List<Callable> callableList = regionTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("Blacklist region analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private class RegionTask implements Callable
    {
        private final ChrBaseRegion mRegion;

        private final List<RepeatMaskData> mRepeatMaskRegions;
        private final List<BaseRegion> mBlackListRegions;
        private final List<UnmappingRegion> mUnmapRegions;

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
        }

        @Override
        public Long call()
        {
            int currentWindowStart = mRegion.start();
            int windowEnd = currentWindowStart + mBaseWindowLength - 1;

            int minRegionEnd = mRegion.end() - (int)(mBaseWindowLength * 0.1);

            while(windowEnd <= minRegionEnd)
            {
                analyseWindow(new BaseRegion(currentWindowStart, min(windowEnd, mRegion.end())));

                currentWindowStart = windowEnd + 1;
                windowEnd = currentWindowStart + mBaseWindowLength - 1;
            }

            return (long)0;
        }

        private void analyseWindow(final BaseRegion region)
        {
            // ignore regions already covered by an unmapping region
            for(UnmappingRegion unmapRegion : mUnmapRegions)
            {
                // check for same required overlap as used in Redux to skip a read
                if(!hasSufficientOverlap(region, unmapRegion, 0.9))
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

            double topBasePerc = topBaseCount / (double)mBaseWindowLength;
            double topAndSecondBasePerc = (topBaseCount + secondBaseCount) / (double)mBaseWindowLength;

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

            writeRegionData(mRegion.Chromosome, region, dominantBases, dominantBasePercent, repeatMasks, blacklistRegions);
        }
    }

    private boolean hasSufficientOverlap(final BaseRegion window, final BaseRegion otherRegion, final double requiredPerc)
    {
        if(!otherRegion.overlaps(window))
            return false;

        int overlap = min(otherRegion.end(), window.end()) - max(otherRegion.start(), window.start()) + 1;
        return overlap >= requiredPerc * mBaseWindowLength;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END);

            sj.add("DominantBases").add("Concentration").add("BlacklistRegions").add("RepeatInfo");

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
            final List<RepeatMaskData> repeatMasks, final List<BaseRegion> blacklistRegions)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(chromosome).add(String.valueOf(region.start())).add(String.valueOf(region.end()));

            sj.add(dominantBases).add(String.format("%.2f", basePercent));

            String blacklistInfo = blacklistRegions.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));

            String repeatInfo = repeatMasks.stream()
                    .map(x -> format("%s:%d-%d", x.ClassType, x.Region.start(), x.Region.end())).collect(Collectors.joining(ITEM_DELIM));

            sj.add(blacklistInfo).add(repeatInfo);

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

        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");
        configBuilder.addInteger(PARTITION_SIZE, "Read length", DEFAULT_CHR_PARTITION_SIZE);

        configBuilder.addInteger(BASE_WINDOW_LENGTH, "Base window length for analysis", DEFAULT_BASE_WINDOW_LENGTH);
        configBuilder.addDecimal(MIN_DOMINANT_PERC, "Min dominant concerntration of 1 or 2 bases", DEFAULT_MIN_DOMINANT_PERC);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BlacklistRepeatAnalyser blacklistRepeatAnalyser = new BlacklistRepeatAnalyser(configBuilder);
        blacklistRepeatAnalyser.run();
    }
}
