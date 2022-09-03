package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.region.GenomeRegion
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator
import htsjdk.samtools.SamReader
import htsjdk.samtools.SamReaderFactory
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException
import java.util.*
import java.util.concurrent.BlockingQueue
import java.util.concurrent.LinkedBlockingDeque

object AsyncBamReader
{
    private const val MAX_BAM_RECORD_Q_SIZE = 100000
    private val logger = LogManager.getLogger(AsyncBamReader::class.java)

    @JvmStatic
    @Throws(InterruptedException::class)
    fun processBam(
        bamFile: String,
        samReaderFactory: SamReaderFactory,
        genomeRegions: Collection<GenomeRegion>,
        asyncRecordHandler: (SAMRecord) -> Unit,
        threadCount: Int
    )
    {
        logger.debug("Processing {} potential sites in bam {}", genomeRegions.size, bamFile)

        // create genome regions from the loci
        // val taskCompletion = BamTaskCompletion(taskQ.size)

        val bamRecordQ: BlockingQueue<Optional<SAMRecord>> = LinkedBlockingDeque()

        // create the bam record consumers
        val recordConsumers = ArrayList<BamRecordConsumerThread>()
        for (i in 0 until Math.max(threadCount, 1))
        {
            val t = BamRecordConsumerThread(bamRecordQ, asyncRecordHandler)
            t.name = String.format("worker-%d", i)
            t.start()
            recordConsumers.add(t)
        }
        logger.info("{} bam record consumer threads started", recordConsumers.size)

        val bamReader = BamReader(bamFile, samReaderFactory, genomeRegions, bamRecordQ)
        bamReader.run()
        bamRecordQ.put(Optional.empty()) // signals consumer to finish

        for (t in recordConsumers)
        {
            while (t.isAlive)
            {
                // check status every 30 seconds
                t.join(30000)
            }
        }

        logger.info("{} bam reader threads finished", recordConsumers.size)
    }

    internal class BamReader(
        bamFile: String,
        samReaderFactory: SamReaderFactory,
        private val genomeRegionList: Collection<GenomeRegion>,
        private val outputRecordQ: Queue<Optional<SAMRecord>>)
    {
        private val mSamReader: SamReader = samReaderFactory.open(File(bamFile))

        fun run()
        {
            logger.debug("bam reader start")
            for (genomeRegion in genomeRegionList)
            {
                mSamReader.queryOverlapping(genomeRegion.chromosome(), genomeRegion.start(), genomeRegion.end())
                    .use({ iterator -> processRecords(iterator) })
            }

            // process unmapped reads
            mSamReader.queryUnmapped().use({ iterator -> processRecords(iterator) })

            try
            {
                mSamReader.close()
            }
            catch (e: IOException)
            {
                logger.error("IO exception in SamReader::close: {}", e.message)
            }
            logger.debug("bam reader finish")
        }

        private fun processRecords(iterator: SAMRecordIterator)
        {
            while (iterator.hasNext())
            {
                val record = iterator.next()

                /*if (    record.getReadName().equals("ST-E00289:125:HLK2HCCXY:6:2205:6512:35098") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:6:2206:30076:31195") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:4:2218:26636:11804") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:7:2207:14539:63384") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:8:2220:8694:54102") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:6:2211:16305:10152") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:6:1202:18243:9888") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:6:1204:20831:45206") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:7:1118:22333:62997") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:7:2119:20892:72174") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:8:1101:15879:59763") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:8:1222:7963:47492") ||
                        record.getReadName().equals("ST-E00289:125:HLK2HCCXY:3:2118:15818:16516"))
                {
                    logger.info("read record is what we are looking for: {}, genome region: {}", record, task.genomeRegion);
                }*/

                if (record.duplicateReadFlag)
                {
                    continue
                }

                // we don't want to check the alignment region
                outputRecordQ.add(Optional.of(record))

                // important: we want to take a break if we are going to flood the output queue
                checkAndWaitForConsumer()
            }
        }

        private fun checkAndWaitForConsumer()
        {
            if (outputRecordQ.size >= MAX_BAM_RECORD_Q_SIZE)
            {
                logger.info("bam record Q size: {}, max size reached, pausing bam reader", outputRecordQ.size)

                // drain at least 80% of the queue before going back to work
                while (outputRecordQ.size * 5 >= MAX_BAM_RECORD_Q_SIZE)
                {
                    try
                    {
                        Thread.sleep(1000)
                    }
                    catch (ex: InterruptedException)
                    {
                        Thread.currentThread().interrupt()
                    }
                }
                logger.info("finished sleeping, bam record Q size: {}", outputRecordQ.size)
            }
        }
    }

    internal class BamRecordConsumerThread(
        private val bamRecordQ: BlockingQueue<Optional<SAMRecord>>,
        private val samRecordHandler: (SAMRecord) -> Unit
    ) : Thread()
    {
        override fun run()
        {
            logger.debug("bam record consumer thread start")
            while (true)
            {
                val record: Optional<SAMRecord> = try
                {
                    bamRecordQ.take()
                }
                catch (e: InterruptedException)
                {
                    break
                }
                if (record.isEmpty)
                {
                    // if record is empty, it signals consumer to stop
                    // we want to put it back so other consumer threads
                    // will also know to stop
                    bamRecordQ.put(record)
                    break
                }
                samRecordHandler(record.get())
            }
            logger.debug("bam record consumer thread finish")
        }
    }
}