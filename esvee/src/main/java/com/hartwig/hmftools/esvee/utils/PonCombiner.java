package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.esvee.caller.annotation.PonCache;
import com.hartwig.hmftools.esvee.caller.annotation.PonSglRegion;
import com.hartwig.hmftools.esvee.caller.annotation.PonSvRegion;

import org.jetbrains.annotations.NotNull;

public class PonCombiner
{
    private final List<PonCache> mPonCaches;
    private final String mOutputSvFile;
    private final String mOutputSglFile;
    private final RefGenomeVersion mRefGenomeVersion;

    private static final String INPUT_SV_PON_FILES = "input_pon_sv_files";
    private static final String INPUT_SGL_PON_FILES = "input_pon_sgl_files";

    private static final String OUTPUT_SV_PON_FILE = "output_pon_sv_file";
    private static final String OUTPUT_SGL_PON_FILE = "output_pon_sgl_file";

    public PonCombiner(final ConfigBuilder configBuilder)
    {
        mPonCaches = Lists.newArrayList();

        String[] svFiles = configBuilder.getValue(INPUT_SV_PON_FILES).split(ITEM_DELIM, -1);
        String[] sglFiles = configBuilder.getValue(INPUT_SGL_PON_FILES).split(ITEM_DELIM, -1);

        if(svFiles.length != sglFiles.length)
        {
            SV_LOGGER.error("requires same number of input files for SVs as for SGLs");
            System.exit(1);
        }

        for(int i = 0; i < svFiles.length; ++i)
        {
            mPonCaches.add(new PonCache(0, 0, svFiles[i], sglFiles[i], true));
        }

        if(mPonCaches.get(0).svRegions().keySet().stream().anyMatch(x -> x.startsWith("chr")))
        {
            mRefGenomeVersion = V38;
        }
        else
        {
            mRefGenomeVersion = V37;
        }

        mOutputSvFile = configBuilder.getValue(OUTPUT_SV_PON_FILE);
        mOutputSglFile = configBuilder.getValue(OUTPUT_SGL_PON_FILE);
    }

    public void run()
    {
        if(mPonCaches.isEmpty() || mOutputSglFile == null || mOutputSvFile == null)
        {
            SV_LOGGER.error("invalid inputs");
            System.exit(1);
        }

        SV_LOGGER.info("SV PON file merge");

        mergeSvPonFiles();
        mergeSglPonFiles();

        SV_LOGGER.info("SV PON merge complete");
    }

    private static final int LOG_COUNT = 1000000;

    private void mergeSvPonFiles()
    {
        SV_LOGGER.info("merging {} SV PON files to: {}", mPonCaches.size(), mOutputSvFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputSvFile, false);

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

                List<PonSvRegion> combinedRegions = Lists.newArrayList();

                for(PonCache ponCache : mPonCaches)
                {
                    List<PonSvRegion> regions = ponCache.svRegions().get(chrStr);

                    if(regions == null)
                        continue;

                    combinedRegions.addAll(regions);
                }

                // sort and merge
                mergeSvRegions(chrStr, combinedRegions);

                SV_LOGGER.debug("chr({}) writing {} SV regions", chrStr, combinedRegions.size());

                for(PonSvRegion region : combinedRegions)
                {
                    writer.write(region.toBedRecord());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write SV PON output file: {}", e.toString());
            System.exit(1);
        }
    }

    @VisibleForTesting
    public static void mergeSvRegions(final String chromosomeStart, final List<PonSvRegion> combinedRegions)
    {
        Collections.sort(combinedRegions);

        if(combinedRegions.size() < 2)
            return;

        int mergeCount = 0;
        int initialCount = combinedRegions.size();

        int index = 0;
        while(index < combinedRegions.size())
        {
            PonSvRegion currentRegion = combinedRegions.get(index);

            while(true)
            {
                boolean foundMerge = false;

                // search all regions starting with this base and merge if both ends overlap
                int nextIndex = index + 1;

                while(nextIndex < combinedRegions.size())
                {
                    PonSvRegion nextRegion = combinedRegions.get(nextIndex);

                    // must overlap at the start
                    if(nextRegion.RegionStart.start() > currentRegion.RegionStart.end())
                        break;

                    if(nextRegion.OrientStart != currentRegion.OrientStart || nextRegion.OrientEnd != currentRegion.OrientEnd)
                    {
                        ++nextIndex;
                        continue;
                    }

                    if(!currentRegion.RegionEnd.overlaps(nextRegion.RegionEnd))
                    {
                        ++nextIndex;
                        continue;
                    }

                    // doesn't matter where the end is - but expand to the longer of the two if any end regions overlap
                    SV_LOGGER.trace("merging region({} -> {}) with next({} -> {})",
                            currentRegion.RegionStart, currentRegion.RegionEnd, nextRegion.RegionStart, nextRegion.RegionEnd);

                    currentRegion.RegionStart.setStart(Math.min(currentRegion.RegionStart.start(), nextRegion.RegionStart.start()));
                    currentRegion.RegionStart.setEnd(max(currentRegion.RegionStart.end(), nextRegion.RegionStart.end()));
                    currentRegion.RegionEnd.setStart(Math.min(currentRegion.RegionEnd.start(), nextRegion.RegionEnd.start()));
                    currentRegion.RegionEnd.setEnd(max(currentRegion.RegionEnd.end(), nextRegion.RegionEnd.end()));

                    currentRegion.setUpdatePonCount(max(currentRegion.PonCount, nextRegion.PonCount));

                    combinedRegions.remove(nextIndex);
                    foundMerge = true;
                    ++mergeCount;
                }

                if(!foundMerge)
                    break;

                // otherwise repeat now region(s) have been expanded
            }

            ++index;
        }

        SV_LOGGER.debug("chr({}) regions initial({}) final({}) merged({})",
                chromosomeStart, initialCount, initialCount - mergeCount, mergeCount);
    }

    private void mergeSglPonFiles()
    {
        SV_LOGGER.info("merging {} SGL PON files to: {}", mPonCaches.size(), mOutputSglFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputSglFile, false);

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

                List<PonSglRegion> combinedRegions = Lists.newArrayList();

                for(PonCache ponCache : mPonCaches)
                {
                    List<PonSglRegion> regions = ponCache.sglRegions().get(chrStr);

                    if(regions == null)
                        continue;

                    combinedRegions.addAll(regions);
                }

                // sort and merge
                mergeSglRegions(chrStr, combinedRegions);

                SV_LOGGER.debug("chr({}) writing {} SGL regions", chrStr, combinedRegions.size());

                for(PonSglRegion ponRegion : combinedRegions)
                {
                    writer.write(ponRegion.toBedRecord());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write SGL PON output file: {}", e.toString());
            System.exit(1);
        }
    }

    @VisibleForTesting
    public static void mergeSglRegions(final String chromosomeStart, final List<PonSglRegion> combinedRegions)
    {
        if(combinedRegions.size() < 2)
            return;

        Collections.sort(combinedRegions);

        int mergeCount = 0;
        int initialCount = combinedRegions.size();

        int index = 0;
        while(index < combinedRegions.size())
        {
            PonSglRegion currentRegion = combinedRegions.get(index);

            // search all regions starting with this base and merge if both ends overlap
            int nextIndex = index + 1;

            while(nextIndex < combinedRegions.size())
            {
                PonSglRegion nextRegion = combinedRegions.get(nextIndex);

                if(nextRegion.Region.start() > currentRegion.Region.end())
                    break;

                if(nextRegion.Orient != currentRegion.Orient)
                {
                    ++nextIndex;
                    continue;
                }

                // doesn't matter where the end is - but expand to the longer of the two if any end regions overlap
                SV_LOGGER.trace("merging region({}) with next({})", currentRegion.Region, nextRegion.Region);

                currentRegion.Region.setEnd(max(currentRegion.Region.end(), nextRegion.Region.end()));
                currentRegion.setUpdatePonCount(max(currentRegion.PonCount, nextRegion.PonCount));

                combinedRegions.remove(nextIndex);
                ++mergeCount;
            }

            ++index;
        }

        SV_LOGGER.debug("chr({}) regions initial({}) final({}) merged({})",
                chromosomeStart, initialCount, initialCount - mergeCount, mergeCount);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(OUTPUT_SGL_PON_FILE, true, "Output SGL PON file");
        configBuilder.addConfigItem(OUTPUT_SV_PON_FILE, true, "Output SV PON file");

        configBuilder.addConfigItem(INPUT_SV_PON_FILES, true, "List of input PON files, separated by ';'");
        configBuilder.addConfigItem(INPUT_SGL_PON_FILES, true, "List of input PON files, separated by ';'");

        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PonCombiner ponCombiner = new PonCombiner(configBuilder);
        ponCombiner.run();
    }
}
