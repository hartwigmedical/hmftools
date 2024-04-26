package com.hartwig.hmftools.bamtools.tofastq;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.logging.log4j.Level;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RemoteReadHandler
{
    private final FastqConfig mConfig;
    private final FastqWriterCache mWriterCache;

    private final HashBamWriter mHashBamWriter;

    private final AtomicInteger mUnmappedReadCount = new AtomicInteger();

    public RemoteReadHandler(final FastqConfig config, final FastqWriterCache writerCache)
    {
        mConfig = config;
        mWriterCache = writerCache;
        mHashBamWriter = new HashBamWriter(config);
    }

    public void cacheRemoteRead(SAMRecord read)
    {
        mHashBamWriter.writeToHashBam(read);
    }

    public void writeRemoteReadsToFastq(ExecutorService executorService, final ThreadData threadData)
            throws ExecutionException, InterruptedException
    {
        // must close and finalise all hashbams writers
        mHashBamWriter.close();

        final List<CompletableFuture<Void>> futures = new ArrayList<>();

        // use multiple threads to process the hash bams
        for(File hashBam : mHashBamWriter.getHashBams())
        {
            futures.add(CompletableFuture.runAsync(() -> processHashBam(threadData, hashBam), executorService));
        }

        // wait for completion
        CompletableFuture.allOf(futures.toArray(CompletableFuture[]::new)).get();

        if(mUnmappedReadCount.get() > 0)
        {
            BT_LOGGER.printf(Level.INFO, "processed %,d total unmatched reads", mUnmappedReadCount.get());
        }
    }

    private void processHashBam(final ThreadData threadData, final File hashBam)
    {
        FastqWriter fastqWriter = threadData.getFastqWriter();

        int numReads = 0;

        Map<String,SAMRecord> unmatchedReads = new HashMap<>();

        try(SamReader samReader = SamReaderFactory.makeDefault().open(hashBam))
        {
            try(final SAMRecordIterator iterator = samReader.iterator())
            {
                while(iterator.hasNext())
                {
                    processSamRecord(iterator.next(), unmatchedReads, fastqWriter);
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
                mWriterCache.writeUnpairedRead(read);
            }
        }

        BT_LOGGER.printf(Level.DEBUG, "hash bam %s, processed %,d unmatched reads", hashBam.getName(), numReads);
    }

    private void processSamRecord(final SAMRecord read, Map<String,SAMRecord> unmatchedReads, FastqWriter fastqWriter)
    {
        if(read.isSecondaryOrSupplementary())
            return;

        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // unexpected
            return;

        if(!read.getReadPairedFlag())
        {
            mWriterCache.processUnpairedRead(read, fastqWriter);
            return;
        }

        mUnmappedReadCount.incrementAndGet();
        SAMRecord mate = unmatchedReads.remove(read.getReadName());

        if(mate != null)
        {
            mWriterCache.processReadPair(read, mate, fastqWriter);
            return;
        }

        unmatchedReads.put(read.getReadName(), read);
    }

    // Write all the unmapped reads into the hash bams
    public void cacheAllUnmappedReads()
    {
        BT_LOGGER.info("start writing unmapped reads to {} hash bams", mHashBamWriter.numHashBams());

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFile)))
        {
            try(final SAMRecordIterator iterator = samReader.queryUnmapped())
            {
                while(iterator.hasNext())
                {
                    final SAMRecord read = iterator.next();

                    if(!read.isSecondaryOrSupplementary() && !read.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                    {
                        mHashBamWriter.writeToHashBam(read);
                    }
                }
            }

            BT_LOGGER.info("finished writing unmapped reads to {} hash bams", mHashBamWriter.numHashBams());
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}
