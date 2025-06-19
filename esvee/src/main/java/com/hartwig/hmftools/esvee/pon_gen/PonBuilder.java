package com.hartwig.hmftools.esvee.pon_gen;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.chrEnd;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.chrStart;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.chromosome;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.locationValues;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.orientEnd;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.orientStart;
import static com.hartwig.hmftools.esvee.pon_gen.PonLocations.orientation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.esvee.caller.annotation.PonSglRegion;
import com.hartwig.hmftools.esvee.caller.annotation.PonSvRegion;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class PonBuilder
{
    private final PonConfig mConfig;
    private final List<String> mSampleIds;

    private final PonStore mPonStore;

    public PonBuilder(final ConfigBuilder configBuilder)
    {
        mConfig = new PonConfig(configBuilder);

        mSampleIds = loadSampleIdsFile(configBuilder);

        mPonStore = new PonStore();
    }

    public void run()
    {
        if(mSampleIds.isEmpty() || mConfig.VcfPath == null)
        {
            SV_LOGGER.error("missing sample IDs file or VCF path in config");
            System.exit(1);
        }

        long startTimeMs = System.currentTimeMillis();

        SV_LOGGER.info("generating PON from {} samples, minSamples({})",
                mSampleIds.size(), mConfig.MinSamples);

        List<String> sampleVcfs = Lists.newArrayListWithCapacity(mSampleIds.size());

        for(String sampleId : mSampleIds)
        {
            String vcfFilename = mConfig.VcfPath.replaceAll("\\*", sampleId);

            if(!Files.exists(Paths.get(vcfFilename)))
            {
                SV_LOGGER.warn("sample({}) missing VCF({})", sampleId, vcfFilename);
                continue;
            }

            sampleVcfs.add(vcfFilename);
        }

        List<PonSampleTask> ponTasks = Lists.newArrayList();

        int threads = min(mConfig.Threads, mSampleIds.size());

        for(int i = 0; i < threads; ++i)
        {
            ponTasks.add(new PonSampleTask(mConfig, mPonStore));
        }

        int taskIndex = 0;
        for(String sampleVcf : sampleVcfs)
        {
            ponTasks.get(taskIndex).addSampleVcf(sampleVcf);
            ++taskIndex;

            if(taskIndex >= threads)
                taskIndex = 0;
        }

        List<Callable> callableList = ponTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, callableList.size()))
            System.exit(1);

        SV_LOGGER.info("building final PON");
        buildSvPON();
        buildSglPON();

        SV_LOGGER.info("Esvee PON building complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private void buildSvPON()
    {
        // convert locations into PON entries, applying filters and merging positions as configured

        Map<String,List<PonSvRegion>> chrPonMap = Maps.newHashMap(); // SV entries keyed by starting/lower chromosome

        for(Map.Entry<String,PonLocations> entry : mPonStore.getSvLocations().entrySet())
        {
            String locationKey = entry.getKey();
            String[] locValues = locationValues(locationKey);

            String chrStart = chrStart(locValues);

            List<PonSvRegion> ponRegions = chrPonMap.get(chrStart);

            if(ponRegions == null)
            {
                ponRegions = Lists.newArrayList();
                chrPonMap.put(chrStart, ponRegions);
            }

            String chrEnd = chrEnd(locValues);
            Orientation orientStart = orientStart(locValues);
            Orientation orientEnd = orientEnd(locValues);

            PonLocations locations = entry.getValue();

            for(LocationCounter locationStart : locations.Locations)
            {
                int positionStart = locationStart.Position;

                ChrBaseRegion regionStart = new ChrBaseRegion(chrStart, positionStart, positionStart);

                for(LocationCounter locationEnd : locationStart.getNextLocations())
                {
                    ChrBaseRegion regionEnd = new ChrBaseRegion(chrEnd, locationEnd.Position, locationEnd.Position);

                    ponRegions.add(new PonSvRegion(regionStart, orientStart, regionEnd, orientEnd, locationEnd.getCount()));
                }
            }
        }

        // sort and combine records
        for(List<PonSvRegion> regions : chrPonMap.values())
        {
            Collections.sort(regions);

            int index = 0;

            while(index < regions.size())
            {
                PonSvRegion region = regions.get(index);

                int nextIndex = index + 1;

                while(nextIndex < regions.size())
                {
                    PonSvRegion nextRegion = regions.get(nextIndex);

                    if(!tryMergeSvRegions(region, nextRegion))
                        break;

                    regions.set(index, mergeSvRegions(region, nextRegion));
                    regions.remove(nextIndex);
                }

                ++index;
            }

            // apply filters
            index = 0;

            while(index < regions.size())
            {
                PonSvRegion region = regions.get(index);

                if(region.PonCount < mConfig.MinSamples)
                    regions.remove(index);
                else
                    ++index;
            }
        }

        writeSvPonFile(chrPonMap);
    }

    private String generateFilename(boolean isSv)
    {
        if(isSv)
            return mConfig.OutputDir + format("sv_pon_%s.%s", mConfig.OutputFilenameSuffix, mConfig.WriteBedFiles ? "bedpe.gz" : "tsv.gz");
        else
            return mConfig.OutputDir + format("sgl_pon_%s.%s", mConfig.OutputFilenameSuffix, mConfig.WriteBedFiles ? "bed.gz" : "tsv.gz");
    }

    private void writeSvPonFile(final Map<String,List<PonSvRegion>> chrPonMap)
    {
        String filename = generateFilename(true);

        SV_LOGGER.info("writing SV PON file: {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            if(!mConfig.WriteBedFiles)
            {
                writer.write(PonSvRegion.header());
                writer.newLine();
            }

            RefGenomeVersion refGenomeVersion = chrPonMap.keySet().stream().anyMatch(x -> x.startsWith(CHR_PREFIX)) ? V38 : V37;

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());

                List<PonSvRegion> regions = chrPonMap.get(chrStr);

                if(regions == null || regions.isEmpty())
                    continue;

                for(PonSvRegion region : regions)
                {
                    if(mConfig.WriteBedFiles)
                        writer.write(region.toBedRecord());
                    else
                        writer.write(region.toTsv());

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write SV PON file: {}", e.toString());
        }
    }

    private boolean tryMergeSvRegions(final PonSvRegion region, final PonSvRegion nextRegion)
    {
        // must match on orientations and have both positions within required distance
        if(!region.RegionEnd.Chromosome.equals(nextRegion.RegionEnd.Chromosome))
            return false;

        if(region.OrientStart != nextRegion.OrientStart || region.OrientEnd != nextRegion.OrientEnd)
            return false;

        int permittedBuffer = mConfig.PositionBuffer;

        if(!areProximateRegions(region.RegionStart, nextRegion.RegionStart, permittedBuffer))
            return false;

        return areProximateRegions(region.RegionEnd, nextRegion.RegionEnd, permittedBuffer);
    }

    private PonSvRegion mergeSvRegions(final PonSvRegion region, final PonSvRegion nextRegion)
    {
        ChrBaseRegion regionStart = new ChrBaseRegion(
                region.RegionStart.Chromosome,
                min(region.RegionStart.start(), nextRegion.RegionStart.start()),
                max(region.RegionStart.end(), nextRegion.RegionStart.end()));

        ChrBaseRegion regionEnd = new ChrBaseRegion(
                region.RegionEnd.Chromosome,
                min(region.RegionEnd.start(), nextRegion.RegionEnd.start()),
                max(region.RegionEnd.end(), nextRegion.RegionEnd.end()));

        return new PonSvRegion(
                regionStart, region.OrientStart, regionEnd, region.OrientEnd, region.PonCount + nextRegion.PonCount);
    }

    private boolean canMergeSglRegions(final PonSglRegion region, final PonSglRegion nextRegion)
    {
        if(region.Orient != nextRegion.Orient)
            return false;

        int permittedBuffer = mConfig.PositionBuffer;

        return areProximateRegions(region.Region, nextRegion.Region, permittedBuffer);
    }

    private PonSglRegion mergeSglRegions(final PonSglRegion region, final PonSglRegion nextRegion)
    {
        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(
                region.Region.Chromosome,
                min(region.Region.start(), nextRegion.Region.start()),
                max(region.Region.end(), nextRegion.Region.end()));

        return new PonSglRegion(chrBaseRegion, region.Orient, region.PonCount + nextRegion.PonCount);
    }

    private static boolean areProximateRegions(final ChrBaseRegion region, final ChrBaseRegion nextRegion, final int permittedBuffer)
    {
        if(positionsOverlap(region.start(), region.end(), nextRegion.start(), nextRegion.end()))
            return true;

        return abs(nextRegion.start() - region.end()) <= permittedBuffer;
    }

    private void buildSglPON()
    {
        Map<String,List<PonSglRegion>> chrPonMap = Maps.newHashMap(); // SGL entries keyed by chromosome

        for(Map.Entry<String,PonLocations> entry : mPonStore.getSglLocations().entrySet())
        {
            String locationKey = entry.getKey();
            String[] locValues = locationValues(locationKey);

            String chromosome = chromosome(locValues);

            List<PonSglRegion> ponRegions = chrPonMap.get(chromosome);

            if(ponRegions == null)
            {
                ponRegions = Lists.newArrayList();
                chrPonMap.put(chromosome, ponRegions);
            }

            Orientation orientation = orientation(locValues);

            PonLocations locations = entry.getValue();

            for(LocationCounter location : locations.Locations)
            {
                ChrBaseRegion region = new ChrBaseRegion(chromosome, location.Position, location.Position);
                ponRegions.add(new PonSglRegion(region, orientation, location.getCount()));
            }
        }

        // sort and combine records
        for(List<PonSglRegion> regions : chrPonMap.values())
        {
            Collections.sort(regions);

            int index = 0;

            while(index < regions.size())
            {
                PonSglRegion region = regions.get(index);

                int nextIndex = index + 1;

                while(nextIndex < regions.size())
                {
                    PonSglRegion nextRegion = regions.get(nextIndex);

                    if(!canMergeSglRegions(region, nextRegion))
                        break;

                    regions.set(index, mergeSglRegions(region, nextRegion));
                    regions.remove(nextIndex);
                }

                ++index;
            }

            // apply filters
            index = 0;

            while(index < regions.size())
            {
                PonSglRegion region = regions.get(index);

                if(region.PonCount < mConfig.MinSamples)
                    regions.remove(index);
                else
                    ++index;
            }
        }

        writeSglPonFile(chrPonMap);
    }

    private void writeSglPonFile(final Map<String,List<PonSglRegion>> chrPonMap)
    {
        String filename = generateFilename(false);

        SV_LOGGER.info("writing SGL PON file: {}", filename);

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            if(!mConfig.WriteBedFiles)
            {
                writer.write(PonSglRegion.header());
                writer.newLine();
            }

            RefGenomeVersion refGenomeVersion = chrPonMap.keySet().stream().anyMatch(x -> x.startsWith(CHR_PREFIX)) ? V38 : V37;

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());

                List<PonSglRegion> regions = chrPonMap.get(chrStr);

                if(regions == null || regions.isEmpty())
                    continue;

                for(PonSglRegion region : regions)
                {
                    if(mConfig.WriteBedFiles)
                        writer.write(region.toBedRecord());
                    else
                        writer.write(region.toTsv());

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write SGL PON file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PonConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        PonBuilder ponBuilder = new PonBuilder(configBuilder);
        ponBuilder.run();
    }
}
