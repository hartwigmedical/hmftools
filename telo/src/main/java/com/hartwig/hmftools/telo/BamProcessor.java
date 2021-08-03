package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;
import static com.hartwig.hmftools.telo.TeloUtils.createPartitions;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingDeque;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;

import htsjdk.samtools.SAMRecord;

// use a blocking queue architecture instead
public class BamProcessor
{
    public static void processBam(TeloConfig config) throws InterruptedException
    {
        final Queue<BaseRegion> baseRegionQ = new ConcurrentLinkedQueue<>();
        final BlockingQueue<TelBamRecord> telBamRecordQ = new LinkedBlockingDeque<>();
        final Set<String> incompleteReadNames = Collections.newSetFromMap(new ConcurrentHashMap<>());
        BamRecordWriter writer = new BamRecordWriter(config, telBamRecordQ, incompleteReadNames);
        writer.setProcessingMateRegions(false);
        final List<BamReader> bamReaders = Lists.newArrayList();

        List<BaseRegion> partitions = createPartitions(config);

        // we use a base region with null chromesome to denote unmapped
        partitions.add(new BaseRegion(null, null));

        // put all into the queue
        baseRegionQ.addAll(partitions);

        // create all the bam readers
        for(int i = 0; i < config.Threads; ++i)
        {
            BamReader bamReader = new BamReader(config, baseRegionQ, telBamRecordQ, incompleteReadNames);
            bamReaders.add(bamReader);
        }

        // start processing threads and run til completion
        runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);

        TE_LOGGER.info("initial BAM file processing complete");

        // now we process the incomplete groups
        List<BaseRegion> mateRegions = getIncompleteReadGroupMateRegions(writer.getmIncompleteReadGroups());

        if(!mateRegions.isEmpty())
        {
            TE_LOGGER.info("processing {} mate regions", mateRegions.size());
            baseRegionQ.addAll(mateRegions);
        }

        writer.setProcessingMateRegions(true);

        // start processing threads and run til completion
        runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);
        writer.writeAllIncompleteReadGroups();
        writer.close();
    }

    private static void runThreadsTillCompletion(Queue<TelBamRecord> telBamRecordQ,
            List<BamReader> bamReaders, BamRecordWriter writer) throws InterruptedException
    {
        List<Thread> bamReaderThreads = bamReaders.stream().map(Thread::new).collect(Collectors.toList());
        bamReaderThreads.forEach(Thread::start);
        Thread writerThread = new Thread(writer);
        writerThread.start();

        // wait for reader to finish
        for(Thread t : bamReaderThreads)
        {
            t.join();
        }

        // now tell writer to finish also
        TelBamRecord telBamRecord = new TelBamRecord();
        telBamRecord.poison = true;
        telBamRecordQ.add(telBamRecord);
        writerThread.join();
    }

    private static List<BaseRegion> getIncompleteReadGroupMateRegions(Map<String, ReadGroup> incompleteReadGroups)
    {
        // TODO - consider organising these into chromosomes and merging any close positions into a single region
        List<BaseRegion> mateRegions = Lists.newArrayList();
        boolean addUnmapped = false;
        for(ReadGroup readGroup : incompleteReadGroups.values())
        {
            for(SAMRecord read : readGroup.Reads)
            {
                if(!read.getReadPairedFlag())
                {
                    continue;
                }

                if(read.getMateUnmappedFlag())
                {
                    addUnmapped = true;
                }
                else
                {
                    mateRegions.add(new BaseRegion(read.getMateReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart() + 1));
                }
            }
        }

        if (addUnmapped)
        {
            mateRegions.add(new BaseRegion(null, null));
        }
        return mateRegions;
    }
}
