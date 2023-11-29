package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.common.region.ChrBaseRegion
import com.hartwig.hmftools.common.samtools.BamSlicer
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealUtils
import com.hartwig.hmftools.teal.TealUtils.createPartitions
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import org.apache.commons.lang3.mutable.MutableInt
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*
import java.util.concurrent.*

// use a blocking queue architecture to process a bam file
// it creates multiple bam reader worker threads to read the bam file, and a writer thread to
// drain the output
object BamProcessor
{
    private val logger = LogManager.getLogger(javaClass)
    private const val REGION_CONSOLIDATION_DISTANCE = 1000
    
    @JvmStatic
    @Throws(InterruptedException::class)
    fun processBam(config: TelbamParams, executorService: ExecutorService)
    {
        // we create one less for the bam writer
        logger.info("processing bam: {}", config.bamFile)
        val telBamRecordQ: BlockingQueue<TelBamRecord> = LinkedBlockingDeque()
        val incompleteReadNames = Collections.newSetFromMap(ConcurrentHashMap<String, Boolean>())
        val writer = BamRecordWriter(config, telBamRecordQ, incompleteReadNames)
        writer.setProcessingMissingReadRegions(false)
        val samReaderList = Collections.synchronizedList(ArrayList<SamReader>())

        // One bam reader per thread instead of one per task
        val threadBamReader = ThreadLocal.withInitial {
            var factory = SamReaderFactory.makeDefault()
            if (config.refGenomeFile != null && !config.refGenomeFile!!.isEmpty())
            {
                factory = factory.referenceSequence(File(config.refGenomeFile!!))
            }
            val samReader: SamReader = factory.open(File(config.bamFile!!))
            samReaderList.add(samReader)
            samReader
        }

        // add unmapped read first cause it is slower to process
        val partitions = ArrayList<ChrBaseRegion?>()
        partitions.add(null)
        partitions.addAll(createPartitions(config))

        // add the bam reader tasks
        val bamReaderTasks: MutableList<Future<*>> = ArrayList()

        for (baseRegion in partitions)
        {
            val task = Runnable { sliceRegionTask(threadBamReader, baseRegion, telBamRecordQ, incompleteReadNames) }
            bamReaderTasks.add(executorService.submit(task))
        }

        // start processing threads and run til completion
        runTasksTillCompletion(telBamRecordQ, bamReaderTasks, writer)
        logger.info("initial BAM file processing complete")

        // we want to process until no new reads have been accepted
        while (!writer.incompleteReadGroups.isEmpty())
        {
            val numAcceptedReads = writer.numAcceptedReads

            // add unmapped read first cause it is slower to process
            partitions.clear()
            partitions.add(null)

            // now we process the incomplete groups, since the threads have already finished there is no
            // concurrency problem
            partitions.addAll(getMissingReadRegions(writer.incompleteReadGroups, config.specificChromosomes))

            logger.info("processing {} missing read regions", partitions.size)

            for (baseRegion in partitions)
            {
                val task = Runnable { sliceRegionTask(threadBamReader, baseRegion, telBamRecordQ, incompleteReadNames) }
                bamReaderTasks.add(executorService.submit(task))
            }
            writer.setProcessingMissingReadRegions(true)

            // start processing threads and run til completion
            runTasksTillCompletion(telBamRecordQ, bamReaderTasks, writer)
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

    private fun sliceRegionTask(
        samReaderSupplier: ThreadLocal<SamReader>,
        region: ChrBaseRegion?,
        telBamRecordQ: Queue<TelBamRecord>,
        incompleteReadNames: Set<String>)
    {
        val reader = samReaderSupplier.get()
        val bamSlicer = BamSlicer(0, true, true, true)
        bamSlicer.setKeepUnmapped()
        val readCount = MutableInt()

        val readHandler = { samRecord: SAMRecord -> processReadRecord(samRecord, telBamRecordQ, incompleteReadNames, readCount) }

        if(region != null)
        {
            //logger.printf(Level.INFO, "processing region(%s:%,d-%,d)", region.chromosome(), region.start(), region.end())
            bamSlicer.slice(reader, region, readHandler)
            logger.printf(Level.INFO, "processed region(%s:%,d-%,d) read count(%s)",
                region.chromosome(), region.start(), region.end(), readCount)
        }
        else
        {
            bamSlicer.queryUnmapped(reader, readHandler)
            logger.info("processed {} unmapped reads", readCount)
        }
    }

    private fun processReadRecord(record: SAMRecord, telBamRecordQ: Queue<TelBamRecord>, incompleteReadNames: Set<String>,
        readCount: MutableInt)
    {
        val hasTeloContent = TealUtils.hasTelomericContent(record.readString)
        if (!hasTeloContent && !incompleteReadNames.contains(record.readName))
        {
            return
        }

        readCount.increment()

        // push it on to the queue
        val telBamRecord = TelBamRecord()
        telBamRecord.samRecord = record
        telBamRecord.hasTeloContent = hasTeloContent
        telBamRecordQ.add(telBamRecord)
    }

    // run the bam reader and record writer threads till completion
    @Throws(InterruptedException::class)
    private fun runTasksTillCompletion(telBamRecordQ: Queue<TelBamRecord>,
                                         bamReaderTasks: MutableList<Future<*>>,
                                         writer: BamRecordWriter)
    {
        val writerThread = Thread(writer)
        writerThread.start()

        // wait for reader tasks to finish
        for (t in bamReaderTasks)
        {
            t.get()
        }
        logger.info("{} bam reader tasks finished", bamReaderTasks.size)

        bamReaderTasks.clear()

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
        return missingReadRegions
    }
}