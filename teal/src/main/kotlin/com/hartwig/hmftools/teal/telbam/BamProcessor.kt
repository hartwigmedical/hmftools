package com.hartwig.hmftools.teal.telbam

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.region.ChrBaseRegion
import com.hartwig.hmftools.common.region.PartitionUtils
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealConstants
import com.hartwig.hmftools.teal.TealUtils.openSamReader
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.samtools.SamReader
import org.apache.logging.log4j.LogManager
import java.util.*
import java.util.concurrent.*

// use a blocking queue architecture to process a bam file
// it creates multiple bam reader worker threads to read the bam file, and a writer thread to
// drain the output
class BamProcessor(val config: TelbamParams)
{
    companion object
    {
        private const val REGION_CONSOLIDATION_DISTANCE = 100_000
    }

    private val logger = LogManager.getLogger(javaClass)

    private val incompleteReadNames: MutableSet<String> = Collections.newSetFromMap(ConcurrentHashMap())
    private val writer = BamRecordWriter(config, incompleteReadNames)
    private val samReaderList: MutableList<SamReader> = Collections.synchronizedList(ArrayList<SamReader>())

    // One bam reader per thread instead of one per task
    private val threadBamReader: ThreadLocal<SamReader> = ThreadLocal.withInitial {
        val samReader: SamReader = openSamReader(config)
        samReaderList.add(samReader)
        samReader
    }

    fun processBam()
    {
        val numBamReaders = Math.max(config.threadCount - 1, 1)
        val numDigits = numBamReaders.toString().length
        val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("thread-%0${numDigits}d").build()
        val executorService = Executors.newFixedThreadPool(numBamReaders, namedThreadFactory)

        // we create one less for the bam writer
        logger.info("processing bam: {}", config.bamFile)
        writer.setProcessingMissingReadRegions(false)

        // add the bam reader tasks
        val bamReaderTasks: MutableList<Future<*>> = ArrayList()

        val partitions = ArrayList<BamPartition>()
        partitions.addAll(createPartitions(config))

        for (partition in partitions)
        {
            val reader = PartitionReader(partition, writer, incompleteReadNames)
            val task = Runnable { reader.readPartition(threadBamReader.get()) }
            bamReaderTasks.add(executorService.submit(task))
        }

        // start processing threads and run til completion
        runTasksTillCompletion(bamReaderTasks)
        logger.info("initial BAM file processing complete")

        // we want to process until no new reads have been accepted
        while (writer.incompleteReadGroups.isNotEmpty())
        {
            val numAcceptedReads = writer.numAcceptedReads

            // add unmapped read first cause it is slower to process
            partitions.clear()

            // now we process the incomplete groups, since the threads have already finished there is no
            // concurrency problem
            partitions.addAll(getMissingReadRegions(writer.incompleteReadGroups, config.specificChromosomes))

            logger.info("processing {} missing read regions", partitions.size)

            for (partition in partitions)
            {
                val reader = PartitionReader(partition, writer, incompleteReadNames)
                val task = Runnable { reader.readPartition(threadBamReader.get()) }
                bamReaderTasks.add(executorService.submit(task))
            }
            writer.setProcessingMissingReadRegions(true)

            // start processing threads and run til completion
            runTasksTillCompletion(bamReaderTasks)
            if (writer.numAcceptedReads == numAcceptedReads)
            {
                // no change, so we did not find any more reads
                break
            }
        }
        writer.finish()

        // close all sam readers
        for (samReader in samReaderList)
        {
            samReader.close()
        }
    }

    // run the bam reader and record writer threads till completion
    @Throws(InterruptedException::class)
    private fun runTasksTillCompletion(bamReaderTasks: MutableList<Future<*>>)
    {
        // wait for reader tasks to finish
        for (t in bamReaderTasks)
        {
            t.get()
        }
        logger.info("{} bam reader tasks finished", bamReaderTasks.size)
        bamReaderTasks.clear()
    }

    private fun getMissingReadRegions(incompleteReadGroups: Map<String, ReadGroup>, specificChromosomes: List<String>): List<BamPartition>
    {
        var missingReadRegions: MutableList<ChrBaseRegion> = ArrayList()
        for (readGroup in incompleteReadGroups.values)
        {
            assert(readGroup.invariant())
            missingReadRegions.addAll(readGroup.findMissingReadBaseRegions())
            if (specificChromosomes.isNotEmpty())
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
        if (missingReadRegions.isNotEmpty())
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
                    if (overlapRegion.end() + REGION_CONSOLIDATION_DISTANCE >= br.start())
                    {
                        // still overlapping, we can set end
                        if(br.end() > overlapRegion.end())
                        {
                            overlapRegion.setEnd(br.end())
                        }
                        continue
                    }
                }
                // no longer can reuse old one, make a new one
                overlapRegion = br.clone()
                compactBaseRegions.add(overlapRegion)
            }
            missingReadRegions = compactBaseRegions
        }

        val partitionList = ArrayList<BamPartition>()
        // add unmapped read first cause it is slower to process
        partitionList.add(BamPartition.ofUnmapped())
        partitionList.addAll(missingReadRegions.map { o -> BamPartition.ofRegion(o) })
        return partitionList
    }

    private fun createPartitions(config: TelbamParams): List<BamPartition>
    {
        val samSequences: List<SAMSequenceRecord>

        openSamReader(config).use { samReader ->
            samSequences = samReader.fileHeader.sequenceDictionary.sequences
        }

        val partitions: MutableList<BamPartition> = ArrayList()

        // add unmapped read first cause it is slower to process
        partitions.add(BamPartition.ofUnmapped())

        val partitionSize = TealConstants.DEFAULT_PARTITION_SIZE
        for (seq in samSequences)
        {
            val chrStr = seq.sequenceName
            if (config.specificChromosomes.isNotEmpty() && !config.specificChromosomes.contains(chrStr)) continue

            partitions.addAll(PartitionUtils.buildPartitions(chrStr, seq.sequenceLength, partitionSize)
                .map { o -> BamPartition.ofRegion(o) })
        }
        return partitions
    }
}