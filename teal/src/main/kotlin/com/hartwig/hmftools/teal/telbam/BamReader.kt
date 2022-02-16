package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.teal.TealUtils.openSamReader
import com.hartwig.hmftools.teal.TealUtils.hasTelomericContent
import java.lang.Runnable
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion
import htsjdk.samtools.SamReader
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import java.util.*

class BamReader(mConfig: TelbamParams,
                private val mTaskQ: Queue<Task>,
                private val mTelBamRecordQ: Queue<TelBamRecord>,
                private val mIncompleteReadNames: Set<String>) : Runnable
{
    private val logger = LogManager.getLogger(javaClass)
    
    class Task
    {
        private var mBaseRegion: ChrBaseRegion? = null
        private var mQueryUnmapped = false
        fun baseRegion(): ChrBaseRegion?
        {
            return mBaseRegion
        }

        fun queryUnmapped(): Boolean
        {
            return mQueryUnmapped
        }

        companion object
        {
            fun fromBaseRegion(baseRegion: ChrBaseRegion?): Task
            {
                val t = Task()
                t.mBaseRegion = baseRegion
                return t
            }

            fun fromQueryUnmapped(): Task
            {
                val t = Task()
                t.mQueryUnmapped = true
                return t
            }
        }
    }

    private var mReadCount = 0

    // one sam reader per thread
    private val mSamReader: SamReader = openSamReader(mConfig)

    override fun run()
    {
        while (true)
        {
            var task: Task
            task = try
            {
                mTaskQ.remove()
            }
            catch (e: NoSuchElementException)
            {
                // finished processing
                break
            }
            findTelomereContent(task)
        }
    }

    fun findTelomereContent(task: Task)
    {
        if (task.queryUnmapped())
        {
            processUnmappedReads()
        }
        if (task.baseRegion() != null)
        {
            processBamByRegion(task.baseRegion())
        }
    }

    // Important note:
    // here we must use SamReader::query. I have tested that this function would return unmapped read that has
    // chromosome and start position set. We need that cause unfortunately the queryUnmapped function is not going to return
    // them. See https://github.com/samtools/htsjdk/issues/278
    private fun processBamByRegion(baseRegion: ChrBaseRegion?)
    {
        logger.info("processing region({})", baseRegion.toString())
        mReadCount = 0
        mSamReader.query(baseRegion!!.Chromosome, baseRegion.start(), baseRegion.end(), false).use { iterator ->
            while (iterator.hasNext())
            {
                val record = iterator.next()
                processReadRecord(record)
            }
        }
        logger.info("processed region({}) read count({})", baseRegion.toString(), mReadCount)
    }

    // Important note:
    // queryUnmapped will not return unmapped read that has chromosome / pos start set. These are
    // unmapped read where the mate is mapped.
    // see https://github.com/samtools/htsjdk/issues/278
    private fun processUnmappedReads()
    {
        logger.info("processing unmapped reads, incomplete={}", mIncompleteReadNames.size)
        mReadCount = 0
        mSamReader.queryUnmapped().use { iterator ->
            while (iterator.hasNext())
            {
                val record = iterator.next()
                processReadRecord(record)
            }
        }
        logger.info("processed {} unmapped reads", mReadCount)
    }

    private fun hasTelomericContent(record: SAMRecord): Boolean
    {
        return hasTelomericContent(record.readString)
    }

    private fun processReadRecord(record: SAMRecord)
    {
        val hasTeloContent = hasTelomericContent(record)
        if (!hasTeloContent && !mIncompleteReadNames.contains(record.readName))
        {
            return
        }
        ++mReadCount

        // push it on to the queue
        val telBamRecord = TelBamRecord()
        telBamRecord.samRecord = record
        telBamRecord.hasTeloContent = hasTeloContent
        mTelBamRecordQ.add(telBamRecord)
    }
}