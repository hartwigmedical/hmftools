package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.teal.TealUtils.createPartitions
import kotlin.Throws
import java.lang.InterruptedException
import java.util.concurrent.ConcurrentLinkedQueue
import java.util.concurrent.BlockingQueue
import java.util.concurrent.ConcurrentHashMap
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion
import com.hartwig.hmftools.teal.ReadGroup
import org.apache.logging.log4j.LogManager
import java.util.*
import java.util.concurrent.LinkedBlockingDeque

// use a blocking queue architecture to process a bam file
// it creates multiple bam reader worker threads to read the bam file, and a writer thread to
// drain the output
object BamProcessor
{
    private val logger = LogManager.getLogger(javaClass)
    
    @JvmStatic
    @Throws(InterruptedException::class)
    fun processBam(config: TelbamParams)
    {
        // we create one less for the bam writer
        val numBamReaders = Math.max(config.threadCount - 1, 1)
        logger.info("processing bam: {} with {} bam reader threads", config.bamFile, numBamReaders)
        val bamReaderTaskQ: Queue<BamReader.Task> = ConcurrentLinkedQueue()
        val telBamRecordQ: BlockingQueue<TelBamRecord> = LinkedBlockingDeque()
        val incompleteReadNames = Collections.newSetFromMap(ConcurrentHashMap<String, Boolean>())
        val writer = BamRecordWriter(config, telBamRecordQ, incompleteReadNames)
        writer.setProcessingMissingReadRegions(false)
        val bamReaders: MutableList<BamReader> = ArrayList()

        // create all the bam readers
        for (i in 0 until numBamReaders)
        {
            val bamReader = BamReader(config, bamReaderTaskQ, telBamRecordQ, Collections.unmodifiableSet(incompleteReadNames))
            bamReaders.add(bamReader)
        }

        // add unmapped read first cause it is slower to process
        bamReaderTaskQ.add(BamReader.Task.fromQueryUnmapped())
        val partitions = createPartitions(config)

        // put all into the queue
        partitions.forEach({ x: ChrBaseRegion? -> bamReaderTaskQ.add(BamReader.Task.fromBaseRegion(x)) })

        // start processing threads and run til completion
        runThreadsTillCompletion(telBamRecordQ, bamReaders, writer)
        logger.info("initial BAM file processing complete")

        // we want to process until no new reads have been accepted
        while (!writer.incompleteReadGroups.isEmpty())
        {
            val numAcceptedReads = writer.numAcceptedReads

            // add unmapped read first cause it is slower to process
            bamReaderTaskQ.add(BamReader.Task.fromQueryUnmapped())

            // now we process the incomplete groups, since the threads have already finished there is no
            // concurrency problem
            val missingReadRegions = getMissingReadRegions(writer.incompleteReadGroups, config.specificChromosomes)
            if (!missingReadRegions.isEmpty())
            {
                logger.info("processing {} missing read regions", missingReadRegions.size)
                missingReadRegions.forEach({ x: ChrBaseRegion? -> bamReaderTaskQ.add(BamReader.Task.fromBaseRegion(x)) })
            }
            writer.setProcessingMissingReadRegions(true)

            // start processing threads and run til completion
            runThreadsTillCompletion(telBamRecordQ, bamReaders, writer)
            if (writer.numAcceptedReads == numAcceptedReads)
            {
                // no change, so we did not find any more reads
                break
            }
        }
        writer.finish()
        check(writer.incompleteReadGroups.isEmpty()) {
            // this should be flagged as an error
            String.format("%d read groups could not be completed", writer.incompleteReadGroups.size)
        }
    }

    // run the bam reader and record writer threads till completion
    @Throws(InterruptedException::class)
    private fun runThreadsTillCompletion(telBamRecordQ: Queue<TelBamRecord>,
                                         bamReaders: List<BamReader>, writer: BamRecordWriter)
    {
        val bamReaderThreads = bamReaders.map({ target: BamReader? -> Thread(target) }).toList()
        bamReaderThreads.forEach({ obj: Thread -> obj.start() })
        val writerThread = Thread(writer)
        writerThread.start()

        // wait for reader to finish
        for (t in bamReaderThreads)
        {
            t.join()
        }
        logger.info("{} bam reader threads finished", bamReaderThreads.size)

        // now tell writer to finish also
        val telBamRecord = TelBamRecord()
        telBamRecord.poison = true
        telBamRecordQ.add(telBamRecord)
        writerThread.join()
        logger.info("writer thread finished")
    }

    private fun getMissingReadRegions(incompleteReadGroups: Map<String, ReadGroup>, specificChromosomes: List<String>): List<ChrBaseRegion>
    {
        var missingReadRegions: MutableList<ChrBaseRegion> = ArrayList()
        for (readGroup in incompleteReadGroups.values)
        {
            assert(readGroup.invariant())
            missingReadRegions.addAll(readGroup.findMissingReadBaseRegions())
            if (!specificChromosomes.isEmpty())
            {
                val list: MutableList<ChrBaseRegion> = ArrayList()
                for (x in missingReadRegions)
                {
                    if (specificChromosomes.contains(x.Chromosome))
                    {
                        list.add(x)
                    }
                }
                missingReadRegions = list
            }
        }

        // sort the mate regions
        missingReadRegions.sort()
        if (!missingReadRegions.isEmpty())
        {
            val compactBaseRegions: MutableList<ChrBaseRegion> = ArrayList()

            // now we go through the mate regions and create a new condense list
            var overlapRegion: ChrBaseRegion? = null
            for (br in missingReadRegions)
            {
                // see if still overlap
                if (overlapRegion != null && br.Chromosome == overlapRegion.Chromosome)
                {
                    // we should be sorted this way
                    assert(br.start() >= overlapRegion.start())
                    if (overlapRegion.end() >= br.start())
                    {
                        // still overlapping, we can set it
                        overlapRegion.setStart(br.start())
                        continue
                    }
                }
                // no longer can reuse old one, make a new one
                overlapRegion = br.clone() as ChrBaseRegion
                compactBaseRegions.add(overlapRegion)
            }
            missingReadRegions = compactBaseRegions
        }
        return missingReadRegions
    }
}