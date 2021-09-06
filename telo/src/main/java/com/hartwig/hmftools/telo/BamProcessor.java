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

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

// use a blocking queue architecture to process a bam file
// it creates multiple bam reader worker threads to read the bam file, and a writer thread to
// drain the output
public class BamProcessor
{
    public static void processBam(TeloConfig config) throws InterruptedException
    {
        TE_LOGGER.info("processing bam: {} with {} threads", config.BamFile, config.Threads);

        final Queue<BamReader.Task> bamReaderTaskQ = new ConcurrentLinkedQueue<>();
        final BlockingQueue<TelBamRecord> telBamRecordQ = new LinkedBlockingDeque<>();
        final Set<String> incompleteReadNames = Collections.newSetFromMap(new ConcurrentHashMap<>());
        BamRecordWriter writer = new BamRecordWriter(config, telBamRecordQ, incompleteReadNames);
        writer.setProcessingMateRegions(false);
        final List<BamReader> bamReaders = new ArrayList<>();

        // create all the bam readers
        // we create one less for the bam writer
        int numBamReaders = Math.max(config.Threads - 1, 1);
        for(int i = 0; i < numBamReaders; ++i)
        {
            BamReader bamReader = new BamReader(config, bamReaderTaskQ, telBamRecordQ, Collections.unmodifiableSet(incompleteReadNames));
            bamReaders.add(bamReader);
        }

        // add unmapped read first cause it is slower to process
        bamReaderTaskQ.add(BamReader.Task.fromQueryUnmapped());

        List<ChrBaseRegion> partitions = createPartitions(config);

        // put all into the queue
        partitions.forEach(x -> bamReaderTaskQ.add(BamReader.Task.fromBaseRegion(x)));

        // start processing threads and run til completion
        runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);

        TE_LOGGER.info("initial BAM file processing complete");

        // we want to process until no new reads have been accepted
        while (!writer.getIncompleteReadGroups().isEmpty())
        {
            int numAcceptedReads = writer.getNumAcceptedReads();

            // add unmapped read first cause it is slower to process
            bamReaderTaskQ.add(BamReader.Task.fromQueryUnmapped());

            // now we process the incomplete groups, since the threads have already finished there is no
            // concurrency problem
            List<ChrBaseRegion> mateRegions = getMissingReadRegions(writer.getIncompleteReadGroups());

            if (!mateRegions.isEmpty())
            {
                TE_LOGGER.info("processing {} mate regions", mateRegions.size());
                mateRegions.forEach(x -> bamReaderTaskQ.add(BamReader.Task.fromBaseRegion(x)));
            }

            writer.setProcessingMateRegions(true);

            // start processing threads and run til completion
            runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);

            if (writer.getNumAcceptedReads() == numAcceptedReads)
            {
                // no change, so we did not find any more reads
                break;
            }
        }

        writer.finish();
    }

    // run the bam reader and record writer threads till completion
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

        TE_LOGGER.info("{} bam reader threads finished", bamReaderThreads.size());

        // now tell writer to finish also
        TelBamRecord telBamRecord = new TelBamRecord();
        telBamRecord.poison = true;
        telBamRecordQ.add(telBamRecord);
        writerThread.join();

        TE_LOGGER.info("writer thread finished");
    }

    private static List<ChrBaseRegion> getMissingReadRegions(Map<String, ReadGroup> incompleteReadGroups)
    {
        List<ChrBaseRegion> missingReadRegions = new ArrayList<>();
        for(ReadGroup readGroup : incompleteReadGroups.values())
        {
            assert(readGroup.invariant());
            missingReadRegions.addAll(readGroup.findMissingReadBaseRegions());
        }

        // sort the mate regions
        Collections.sort(missingReadRegions);

        if (!missingReadRegions.isEmpty())
        {
            List<ChrBaseRegion> compactMateRegions = new ArrayList<>();

            // now we go through the mate regions and create a new condense list
            ChrBaseRegion overlapRegion = null;
            for(ChrBaseRegion br : missingReadRegions)
            {
                // see if still overlap
                if(overlapRegion != null && br.Chromosome.equals(overlapRegion.Chromosome))
                {
                    // we should be sorted this way
                    assert(br.start() >= overlapRegion.start());
                    if(overlapRegion.end() >= br.start())
                    {
                        // still overlapping, we can set it
                        overlapRegion.setStart(br.start());
                        continue;
                    }
                }
                // no longer can reuse old one, make a new one
                overlapRegion = (ChrBaseRegion)br.clone();
                compactMateRegions.add(overlapRegion);
            }

            missingReadRegions = compactMateRegions;
        }

        return missingReadRegions;
    }
}
