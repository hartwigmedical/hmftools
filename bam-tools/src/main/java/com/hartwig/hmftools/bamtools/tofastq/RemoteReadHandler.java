package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.atomic.AtomicLong;

import com.hartwig.hmftools.bamtools.common.HashBamWriter;

import org.apache.logging.log4j.Level;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RemoteReadHandler
{
    private static final int UNMAPPED_READ_CHUNK_SIZE = 100_000;
    private static final int NUM_HASH_BAMS = 256;

    private final ToFastqConfig mConfig;

    private final HashBamWriter mHashBamWriter;

    private final AtomicLong mRemoteReadCount = new AtomicLong();

    public RemoteReadHandler(final ToFastqConfig config)
    {
        mConfig = config;

        SAMFileHeader samHeader;
        try(SamReader samReader = ToFastqUtils.openSamReader(mConfig))
        {
            samHeader = samReader.getFileHeader();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        String tempDirPrefix = Paths.get(mConfig.BamFile).getFileName().toString().replace('.', '_') + "_hashbams_";
        mHashBamWriter = new HashBamWriter(samHeader, tempDirPrefix, NUM_HASH_BAMS);
    }

    public void cacheRemoteRead(SAMRecord read)
    {
        mHashBamWriter.writeToHashBam(read);
    }

    public void writeRemoteReadsToFastq(ExecutorService executorService, final ThreadData threadData)
            throws ExecutionException, InterruptedException
    {
        // must close and finalise all hashbam writers
        mHashBamWriter.close();

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        // use multiple threads to process the hash bams
        for(File hashBam : mHashBamWriter.getHashBams().values())
        {
            futures.add(CompletableFuture.runAsync(() -> processHashBam(threadData, hashBam), executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        BT_LOGGER.printf(Level.INFO, "processed %,d total remote reads", mRemoteReadCount.get());
    }

    private void processHashBam(final ThreadData threadData, final File hashBam)
    {
        FastqWriterCache fastqWriterCache = threadData.getFastqWriterCache();

        int numReads = 0;

        Map<String,SAMRecord> unmatchedReads = new HashMap<>();

        try(SamReader samReader = SamReaderFactory.makeDefault().open(hashBam))
        {
            try(final SAMRecordIterator iterator = samReader.iterator())
            {
                while(iterator.hasNext())
                {
                    processSamRecord(iterator.next(), unmatchedReads, fastqWriterCache);
                    ++numReads;
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }

        // now process all remaining reads that are not paired. Log error if they are paired but mate
        // cannot be found
        if(!unmatchedReads.isEmpty())
        {
            BT_LOGGER.info("writing {} unmapped & unpaired reads", unmatchedReads.size());
            for(SAMRecord read : unmatchedReads.values())
            {
                if(read.getReadPairedFlag())
                {
                    BT_LOGGER.error("mate not found for paired read: {}", read);
                }
                fastqWriterCache.writeUnpairedRead(read);
            }
        }

        BT_LOGGER.printf(Level.DEBUG, "processed %,d remote reads in %s", numReads, hashBam.getName());
    }

    private void processSamRecord(final SAMRecord read, Map<String,SAMRecord> unmatchedReads, FastqWriterCache fastqWriterCache)
    {
        if(ToFastqUtils.canIgnoreRead(read))
            return;

        // check for hard clip
        if(read.getCigar().containsOperator(CigarOperator.HARD_CLIP))
        {
            BT_LOGGER.error("read: {}, hard clip found, require extra logic to handle", read);
            throw new RuntimeException("hard clip found on read");
        }

        if(!read.getReadPairedFlag())
        {
            fastqWriterCache.writeUnpairedRead(read);
            return;
        }

        mRemoteReadCount.incrementAndGet();
        SAMRecord mate = unmatchedReads.remove(read.getReadName());

        if(mate != null)
        {
            fastqWriterCache.writeReadPair(read, mate);
            return;
        }

        unmatchedReads.put(read.getReadName(), read);
    }

    // Write all the unmapped reads into the hash bams
    public void cacheAllUnmappedReads(int numTasks, int taskId)
    {
        BT_LOGGER.info("start writing unmapped reads to {} hash bams (task {} of {})",
                mHashBamWriter.numHashBams(), taskId, numTasks);

        try(SamReader samReader = ToFastqUtils.openSamReader(mConfig))
        {
            int readCount = 0;
            try(final SAMRecordIterator iterator = samReader.queryUnmapped())
            {
                while(iterator.hasNext())
                {
                    final SAMRecord read = iterator.next();
                    readCount++;

                    if((readCount / UNMAPPED_READ_CHUNK_SIZE) % numTasks != taskId)
                    {
                        // we divide the reads by 100k chunks between threads
                        continue;
                    }

                    if(!read.isSecondaryOrSupplementary() && !read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                    {
                        mHashBamWriter.writeToHashBam(read);
                    }
                }
            }

            BT_LOGGER.printf(Level.INFO, "finished writing %,d unmapped reads to %d hash bams (task %d of %d)",
                    readCount, mHashBamWriter.numHashBams(), taskId, numTasks);
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}
