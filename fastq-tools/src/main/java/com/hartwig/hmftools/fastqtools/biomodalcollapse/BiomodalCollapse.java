package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static java.lang.String.format;

import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalCollapseUtil.nextFastqRecord;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import com.google.common.collect.Lists;

import htsjdk.samtools.fastq.FastqRecord;

public class BiomodalCollapse
{
    private final BiomodalCollapseConfig mConfig;

    public BiomodalCollapse(final BiomodalCollapseConfig config)
    {
        mConfig = config;
    }

    private Map<String, FastqRecord> loadRefResolvedFastq()
    {
        if(mConfig.RefResolvedFastqPath == null)
        {
            return Collections.emptyMap();
        }

        Map<String, FastqRecord> refResolvesFastqMap = Maps.newHashMap();
        try(BufferedReader fastqReader = createBufferedReader(mConfig.RefResolvedFastqPath))
        {
            FastqRecord fastq = nextFastqRecord(fastqReader);
            while(fastq != null)
            {
                refResolvesFastqMap.put(fastq.getReadName(), fastq);
                fastq = nextFastqRecord(fastqReader);
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }

        return refResolvesFastqMap;
    }

    public void run()
    {
        FQ_LOGGER.info("starting BimodalCollapse");
        long startTimeMs = System.currentTimeMillis();

        Map<String, FastqRecord> refResolvesFastqMap = loadRefResolvedFastq();

        BiomodalCollapseStats stats = new BiomodalCollapseStats();
        List<BiomodalCollapseWorker> workers = Lists.newArrayList();

        BufferedReader fastq1Reader = null;
        BufferedReader fastq2Reader = null;
        BufferedWriter resolvedFastqWriter = null;
        BufferedWriter debugStatsWriter = null;
        try
        {
            fastq1Reader = createBufferedReader(mConfig.Fastq1Path);
            fastq2Reader = createBufferedReader(mConfig.Fastq2Path);
            resolvedFastqWriter = createBufferedWriter(mConfig.CollapsedFastqOutputPath);
            debugStatsWriter = mConfig.DebugStatsOutputPath == null ? null : createBufferedWriter(mConfig.DebugStatsOutputPath);

            if(debugStatsWriter != null)
            {
                debugStatsWriter.write(Arrays.stream(BiomodalConstants.STAT_HEADERS).collect(Collectors.joining(BiomodalConstants.STAT_DELIMITER)));
                debugStatsWriter.newLine();
            }

            SynchronizedPairedFastqReader fastqPairReader =
                    new SynchronizedPairedFastqReader(fastq1Reader, fastq2Reader, mConfig.MaxFastqPairsProcessed);

            if(mConfig.Threads == 1)
            {
                BiomodalCollapseWorker worker =
                        new BiomodalCollapseWorker(0, fastqPairReader, resolvedFastqWriter, debugStatsWriter, refResolvesFastqMap, stats);
                worker.run();
                workers.add(worker);
            }
            else
            {
                for(int i = 0; i < mConfig.Threads; i++)
                {
                    BiomodalCollapseWorker worker =
                            new BiomodalCollapseWorker(i, fastqPairReader, resolvedFastqWriter, debugStatsWriter, refResolvesFastqMap, stats);
                    worker.start();
                    workers.add(worker);
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
        finally
        {

            if(workers.size() > 1)
            {
                for(BiomodalCollapseWorker worker : workers)
                {
                    try
                    {
                        worker.join();
                    }
                    catch(InterruptedException e)
                    {
                        FQ_LOGGER.warn("Failed to join worker thread: threadId({}) {}", worker.ThreadId, e);
                    }
                }
            }

            closeBufferedReader(fastq1Reader);
            closeBufferedReader(fastq2Reader);
            closeBufferedWriter(resolvedFastqWriter);
            closeBufferedWriter(debugStatsWriter);
        }

        FQ_LOGGER.info(stats.toString());
        FQ_LOGGER.info("BimodalCollapse complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public synchronized static void writeResolvedFastqRecord(final BufferedWriter writer, final FastqRecord resolvedFastq)
    {
        StringBuilder readStr = new StringBuilder(resolvedFastq.getReadString());
        StringBuilder MMTagSkipsStr = new StringBuilder();
        int skip = 0;
        for(int i = 0; i < readStr.length(); i++)
        {
            char base = readStr.charAt(i);
            if(base == 'C')
            {
                skip++;
            }
            else if(base == (char) BiomodalConstants.MODC_BASE)
            {
                readStr.setCharAt(i, 'C');
                MMTagSkipsStr.append(',');
                MMTagSkipsStr.append(skip);
                skip = 0;
            }
        }

        try
        {
            writer.write('@');
            writer.write(resolvedFastq.getReadName());
            writer.write('\t');
            writer.write(format("MM:Z:C+C.%s;", MMTagSkipsStr));
            writer.newLine();
            writer.write(readStr.toString());
            writer.newLine();
            writer.write(resolvedFastq.getBaseQualityHeader());
            writer.newLine();
            writer.write(resolvedFastq.getBaseQualityString());
            writer.newLine();
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public synchronized static void writeStatLine(final BufferedWriter writer, final String statLine)
    {
        try
        {
            writer.write(statLine);
            writer.newLine();
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        BiomodalCollapseConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        // set default exception handler for all threads
        // must do this otherwise unhandled exception in other threads might not be reported
        Thread.setDefaultUncaughtExceptionHandler((Thread t, Throwable e) ->
        {
            FQ_LOGGER.fatal("[{}]: uncaught exception: {}", t, e);
            e.printStackTrace(System.err);
            System.exit(1);
        });

        BiomodalCollapseConfig config = new BiomodalCollapseConfig(configBuilder);
        BiomodalCollapse bimodalCollapse = new BiomodalCollapse(config);
        bimodalCollapse.run();
    }
}
