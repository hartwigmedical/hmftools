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

import htsjdk.samtools.SAMRecord;

// use a blocking queue architecture to process a bam file
// it creates multiple bam reader worker threads to read the bam file, and a writer thread to
// drain the output
public class BamProcessor
{
    public static void processBam(TeloConfig config) throws InterruptedException
    {
        final Queue<ChrBaseRegion> baseRegionQ = new ConcurrentLinkedQueue<>();
        final BlockingQueue<TelBamRecord> telBamRecordQ = new LinkedBlockingDeque<>();
        final Set<String> incompleteReadNames = Collections.newSetFromMap(new ConcurrentHashMap<>());
        BamRecordWriter writer = new BamRecordWriter(config, telBamRecordQ, incompleteReadNames);
        writer.setProcessingMateRegions(false);
        final List<BamReader> bamReaders = new ArrayList<>();

        List<ChrBaseRegion> partitions = createPartitions(config);

        // add unmapped base region, add it first cause it is slower to process
        partitions.add(0, TeloConstants.UNMAPPED_BASE_REGION);

        // put all into the queue
        baseRegionQ.addAll(partitions);

        // create all the bam readers
        for(int i = 0; i < config.Threads; ++i)
        {
            BamReader bamReader = new BamReader(config, baseRegionQ, telBamRecordQ, Collections.unmodifiableSet(incompleteReadNames));
            bamReaders.add(bamReader);
        }

        // start processing threads and run til completion
        runThreadsTillCompletion(telBamRecordQ, bamReaders, writer);

        TE_LOGGER.info("initial BAM file processing complete");

        // now we process the incomplete groups, sicne the threads have already finished there is no
        // concurrency problem
        List<ChrBaseRegion> mateRegions = getIncompleteReadGroupMateRegions(writer.getmIncompleteReadGroups());

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

    private static List<ChrBaseRegion> getIncompleteReadGroupMateRegions(Map<String, ReadGroup> incompleteReadGroups)
    {
        List<ChrBaseRegion> mateRegions = new ArrayList<>();
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
                    mateRegions.add(new ChrBaseRegion(read.getMateReferenceName(), read.getMateAlignmentStart(), read.getMateAlignmentStart()));
                }
            }
        }

        // sort the mate regions
        Collections.sort(mateRegions);

        if (!mateRegions.isEmpty())
        {
            List<ChrBaseRegion> compactMateRegions = new ArrayList<>();

            // now we go through the mate regions and create a new condense list
            ChrBaseRegion overlapRegion = null;
            for(ChrBaseRegion br : mateRegions)
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

            mateRegions = compactMateRegions;
        }

        if (addUnmapped)
        {
            mateRegions.add(0, TeloConstants.UNMAPPED_BASE_REGION);
        }
        return mateRegions;
    }
}
