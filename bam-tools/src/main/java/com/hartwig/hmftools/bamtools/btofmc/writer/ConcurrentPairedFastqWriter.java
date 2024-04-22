package com.hartwig.hmftools.bamtools.btofmc.writer;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig.BFQ_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;

import com.google.common.collect.Lists;
import com.google.common.io.CharSink;
import com.google.common.io.Closer;
import com.google.common.io.Files;
import com.google.inject.Inject;
import com.hartwig.hmftools.bamtools.btofmc.BamToFastqConfig;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.apache.commons.lang3.tuple.Pair;
import org.pcollections.PSortedMap;
import org.pcollections.TreePMap;

import htsjdk.samtools.SAMRecord;

// TODO NEXT: TEST
// TODO: Mock file system.
public class ConcurrentPairedFastqWriter extends PairedFastqWriterInterface
{
    private final BamToFastqConfig mConfig;
    // TODO: Wrap into AtomicPMap?
    private final AtomicReference<PSortedMap<Long, NoSyncPairedFastqWriter>> mWritersByThreadRef;
    private final Closer mCloser;

    private final PerformanceCounter mMergePc;

    @Inject
    public ConcurrentPairedFastqWriter(final BamToFastqConfig config)
    {
        mConfig = config;
        mWritersByThreadRef = new AtomicReference<>(TreePMap.empty());
        mCloser = Closer.create();

        mMergePc = new PerformanceCounter("PairedFastqMerge");
    }

    private Pair<String, String> threadFastqFilenames(long threadId)
    {
        String fastq1Filename = format("%s.part_%d.tmp", mConfig.Fastq1OutputFile, threadId);
        String fastq2Filename = format("%s.part_%d.tmp", mConfig.Fastq2OutputFile, threadId);
        return Pair.of(fastq1Filename, fastq2Filename);
    }

    private NoSyncPairedFastqWriter createThreadPairedFastqWriter()
    {
        if(mConfig.NoWrite)
        {
            return mCloser.register(new NoSyncPairedFastqWriter(mConfig, null, null));
        }

        Pair<String, String> fastqFilenames = threadFastqFilenames(Thread.currentThread().getId());
        CharSink fastq1 = Files.asCharSink(new File(fastqFilenames.getLeft()), StandardCharsets.UTF_8);
        CharSink fastq2 = Files.asCharSink(new File(fastqFilenames.getRight()), StandardCharsets.UTF_8);
        return mCloser.register(new NoSyncPairedFastqWriter(mConfig, fastq1, fastq2));
    }

    private NoSyncPairedFastqWriter threadWriter()
    {
        long threadId = Thread.currentThread().getId();
        NoSyncPairedFastqWriter writer = mWritersByThreadRef.get().get(threadId);
        if(writer != null)
        {
            return writer;
        }

        return mWritersByThreadRef.updateAndGet(ref ->
        {
            if(ref.containsKey(threadId))
            {
                return ref;
            }

            return ref.plus(threadId, createThreadPairedFastqWriter());
        }).get(threadId);
    }

    @Override
    public void writePairedFastqRecord(final SAMRecord read1, final SAMRecord read2)
    {
        PairedFastqWriterInterface writer = threadWriter();
        writer.writePairedFastqRecord(read1, read2);
    }

    private void mergeFiles()
    {
        BFQ_LOGGER.info("Merging per thread fastq files");

        // TODO: Move away from cat.
        //        for(long threadId : mWritersByThreadRef.get().keySet())
        //        {
        //            Pair<String, String> fastqFilenames = threadFastqFilenames(threadId);
        //            try(
        //                    BufferedReader fastq1Reader = createBufferedReader(fastqFilenames.getLeft());
        //                    BufferedReader fastq2Reader = createBufferedReader(fastqFilenames.getRight()))
        //            {
        //                while (true)
        //                {
        //                    String line1 = fastq1Reader.readLine();
        //                    String line2 = fastq2Reader.readLine();
        //
        //                    if (line1 == null && line2 == null)
        //                    {
        //                        break;
        //                    }
        //
        //                    if (line1 == null || line2 == null)
        //                    {
        //                        BFQ_LOGGER.error("{} and {} have a differing number of lines", fastqFilenames.getLeft(), fastqFilenames.getRight());
        //                        System.exit(1);
        //                    }
        //
        //                    mFastq1Writer.write(line1);
        //                    mFastq1Writer.newLine();
        //
        //                    mFastq2Writer.write(line2);
        //                    mFastq2Writer.newLine();
        //                }
        //            }
        //        }

        List<String> fastq1Files = Lists.newArrayList();
        List<String> fastq2Files = Lists.newArrayList();
        for(long threadId : mWritersByThreadRef.get().keySet())
        {
            Pair<String, String> fastqFilenames = threadFastqFilenames(threadId);
            fastq1Files.add(fastqFilenames.getLeft());
            fastq2Files.add(fastqFilenames.getRight());
        }

        try
        {
            Process p1 = catMerge(fastq1Files, mConfig.Fastq1OutputFile);
            Process p2 = catMerge(fastq2Files, mConfig.Fastq2OutputFile);
            p1.waitFor();
            p2.waitFor();
        }
        catch(InterruptedException e)
        {
            throw new RuntimeException("Cannot merge per thread fastq output files", e);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Cannot merge per thread fastq output files", e);
        }
    }

    // TODO: Move way from cat merge.
    private static Process catMerge(final List<String> inFiles, final String outfile) throws IOException
    {
        // TODO: Hardcoded path
        final String catPath = "/bin/cat";

        final List<String> commands = Lists.newArrayList();
        commands.add(catPath);
        commands.addAll(inFiles);

        ProcessBuilder processBuilder = new ProcessBuilder(commands);
        // TODO: check that this doesn't append
        processBuilder.redirectOutput(new File(outfile));

        return processBuilder.start();
    }

    private void deleteTempFiles()
    {
        BFQ_LOGGER.info("Deleting temporary per thread fastq files");
        for(long threadId : mWritersByThreadRef.get().keySet())
        {
            Pair<String, String> fastqFilenames = threadFastqFilenames(threadId);
            File tmpFile = new File(fastqFilenames.getLeft());
            tmpFile.delete();
            tmpFile = new File(fastqFilenames.getRight());
            tmpFile.delete();
        }
    }

    @Override
    public void close() throws IOException
    {
        mCloser.close();

        if(mConfig.NoWrite)
        {
            return;
        }

        mMergePc.start();
        mergeFiles();
        mMergePc.stop();

        deleteTempFiles();
    }

    @Override
    public void logStats()
    {
        NoSyncPairedFastqWriter firstThreadWriter = null;
        for(NoSyncPairedFastqWriter threadWriter : mWritersByThreadRef.get().values())
        {
            if(firstThreadWriter == null)
            {
                firstThreadWriter = threadWriter;
                continue;
            }

            firstThreadWriter.mergeStats(threadWriter);
        }

        firstThreadWriter.logStats();
        if(mConfig.PerfDebug)
        {
            mMergePc.logStats();
        }
    }
}
