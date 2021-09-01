package com.hartwig.hmftools.svtools.pon;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.pon.PonLocations.chrEnd;
import static com.hartwig.hmftools.svtools.pon.PonLocations.chrStart;
import static com.hartwig.hmftools.svtools.pon.PonLocations.chromosome;
import static com.hartwig.hmftools.svtools.pon.PonLocations.locationValues;
import static com.hartwig.hmftools.svtools.pon.PonLocations.orientEnd;
import static com.hartwig.hmftools.svtools.pon.PonLocations.orientStart;
import static com.hartwig.hmftools.svtools.pon.PonLocations.orientation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PonBuilder
{
    private final PonConfig mConfig;
    private final Map<String,String> mSampleVcfFiles;
    private final PonStore mPonStore;

    public static final Logger PON_LOGGER = LogManager.getLogger(PonBuilder.class);

    public PonBuilder(final CommandLine cmd)
    {
        mConfig = new PonConfig(cmd);

        mSampleVcfFiles = Maps.newHashMap();
        mPonStore = new PonStore();
    }

    public void run()
    {
        PON_LOGGER.info("building PON from {} samples", mConfig.SampleIds.size());

        findVcfFiles();

        if(mSampleVcfFiles.isEmpty())
        {
            PON_LOGGER.error("missing VCF or batch-run directory");
            System.exit(1);
        }

        PON_LOGGER.info("processing {} samples", mSampleVcfFiles.size());

        List<PonSampleTask> ponTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            int threads = min(mConfig.Threads, mConfig.SampleIds.size());

            for(int i = 0; i < threads; ++i)
            {
                ponTasks.add(new PonSampleTask(i, mPonStore));
            }

            int taskIndex = 0;
            for (final String sampleId : mConfig.SampleIds)
            {
                ponTasks.get(taskIndex).getSampleVcfFiles().put(sampleId, mSampleVcfFiles.get(sampleId));
                ++taskIndex;

                if(taskIndex >= threads)
                    taskIndex = 0;
            }

            final List<Callable> callableList = ponTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            PonSampleTask ponTask = new PonSampleTask(0, mPonStore);
            ponTask.getSampleVcfFiles().putAll(mSampleVcfFiles);
            ponTask.processSamples();
        }

        PON_LOGGER.info("writing PON files");

        writePonFiles();

        PON_LOGGER.info("PON creation complete");
    }

    private void findVcfFiles()
    {
        if(!mConfig.SampleVcfFiles.isEmpty())
        {
            mSampleVcfFiles.putAll(mConfig.SampleVcfFiles);
            return;
        }

        for(String sampleId : mConfig.SampleIds)
        {
            for(String vcfPathPattern : mConfig.VcfFilePatterns)
            {
                String sampleVcfFile = vcfPathPattern.replaceAll("\\*", sampleId);

                if(Files.exists(Paths.get(sampleVcfFile)))
                {
                    mSampleVcfFiles.put(sampleId, sampleVcfFile);
                    break;
                }
            }

            if(!mSampleVcfFiles.containsKey(sampleId))
            {
                PON_LOGGER.warn("sample({}) no VCF file found", sampleId);
            }
        }

        PON_LOGGER.info("found {} VCF files out of {} samples", mSampleVcfFiles.size(), mConfig.SampleIds.size());
    }

    private void writePonFiles()
    {
        try
        {
            String svPonFile = mConfig.OutputDir + "sv_pon.csv";
            BufferedWriter svWriter = createBufferedWriter(svPonFile, false);
            svWriter.write("ChrStart,OrientStart,PosStart,ChrEnd,OrientEnd,PosEnd,PonCount");
            svWriter.newLine();

            for(Map.Entry<String,PonLocations> entry : mPonStore.getSvLocations().entrySet())
            {
                String locationKey = entry.getKey();
                String[] locValues = locationValues(locationKey);
                PonLocations locations = entry.getValue();

                for(LocationCounter locationStart : locations.Locations)
                {
                    for(LocationCounter locationEnd : locationStart.getNextLocations())
                    {
                        if(locationEnd.getCount() < mConfig.MinPonWriteCount)
                            continue;

                        svWriter.write(String.format("%s,%d,%d,%s,%d,%d,%d",
                                chrStart(locValues), orientStart(locValues), locationStart.Position,
                                chrEnd(locValues), orientEnd(locValues), locationEnd.Position, locationEnd.getCount()));
                        svWriter.newLine();
                    }
                }
            }

            svWriter.close();

            String sglPonFile = mConfig.OutputDir + "sgl_pon.csv";
            BufferedWriter sglWriter = createBufferedWriter(sglPonFile, false);
            sglWriter.write("Chromosome,Orientation,Position,PonCount");
            sglWriter.newLine();

            for(Map.Entry<String,PonLocations> entry : mPonStore.getSglLocations().entrySet())
            {
                String locationKey = entry.getKey();
                String[] locValues = locationValues(locationKey);
                PonLocations locations = entry.getValue();

                for(LocationCounter locationStart : locations.Locations)
                {
                    if(locationStart.getCount() < mConfig.MinPonWriteCount)
                        continue;

                    sglWriter.write(String.format("%s,%d,%d,%d",
                            chromosome(locValues), orientation(locValues), locationStart.Position, locationStart.getCount()));
                    sglWriter.newLine();
                }
            }

            sglWriter.close();
        }
        catch(IOException e)
        {
            PON_LOGGER.error("failed to write PON files: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        PonConfig.addOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PonBuilder germlineVcfReader = new PonBuilder(cmd);
        germlineVcfReader.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
