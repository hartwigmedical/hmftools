package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.common.bam.SamRecordUtils
import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealUtils.hasTelomericContent
import com.hartwig.hmftools.teal.TealUtils.openSamReader
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMTag
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import tech.tablesaw.api.BooleanColumn
import tech.tablesaw.api.IntColumn
import tech.tablesaw.api.StringColumn
import tech.tablesaw.api.Table
import tech.tablesaw.io.csv.CsvWriteOptions
import java.util.zip.GZIPOutputStream

class BamRecordWriter(config: TelbamParams, private val mIncompleteReadNames: MutableSet<String>)
{
    private val logger = LogManager.getLogger(javaClass)
    
    private val mIncompleteReadGroups: MutableMap<String, ReadGroup> = HashMap()
    private val mReadDataTable = Table.create("ReadData")
    private val mBamFileWriter: SAMFileWriter
    private var mNumCompletedGroups = 0
    var numAcceptedReads = 0
        private set

    @Volatile
    private var mProcessingMateRegions = false
    private val mReadDataCsvPath: String?

    val incompleteReadGroups: Map<String, ReadGroup>
        get() = mIncompleteReadGroups

    @Synchronized fun setProcessingMissingReadRegions(b: Boolean)
    {
        mProcessingMateRegions = b
    }

    @Synchronized fun finish()
    {
        // we write out the final incomplete group
        mIncompleteReadGroups.values.forEach({ rg -> writeReadGroup(rg) })
        for (rg in mIncompleteReadGroups.values)
        {
            logger.warn("incomplete read group: id={}, is complete={}", rg.name, rg.isComplete(Level.WARN))
            for (r in rg.Reads)
            {
                logger.warn(
                    "record({}) cigar({}) neg strand({}) suppl flag({}) suppl attr({})",
                    r, r.cigarString, r.readNegativeStrandFlag, r.supplementaryAlignmentFlag,
                    r.getStringAttribute(SAMTag.SA.name))
            }
            for (r in rg.SupplementaryReads)
            {
                logger.warn(
                    "supplementary({}) cigar({}) neg strand({}) suppl flag({}) suppl attr({})",
                    r, r.cigarString, r.readNegativeStrandFlag, r.supplementaryAlignmentFlag,
                    r.getStringAttribute(SAMTag.SA.name))
            }
        }
        mBamFileWriter.close()
        writeReadDataToTsv()
        logger.info(
            "wrote {} read groups, complete({}), incomplete({})",
            mNumCompletedGroups + mIncompleteReadGroups.size,
            mNumCompletedGroups, mIncompleteReadGroups.size)
    }

    @Synchronized fun processReadRecord(record: SAMRecord, hasTelomereContent: Boolean)
    {
        var readGroup = mIncompleteReadGroups[record.readName]
        if (readGroup == null)
        {
            if (mProcessingMateRegions)
            {
                // do not create a new group if we are doing mate regions
                // the group would have already been created and written out
                return
            }
            if (hasTelomereContent)
            {
                // cache if new
                readGroup = ReadGroup(record.readName)
                mIncompleteReadGroups[record.readName] = readGroup
                mIncompleteReadNames.add(readGroup.name)
            }
            else
            {
                return
            }
        }
        if (!readGroup.contains(record))
        {
            acceptRead(readGroup, record)
        }
        assert(readGroup.invariant())
        if (readGroup.isComplete()) processCompleteReadGroup(readGroup)
    }

    private fun acceptRead(readGroup: ReadGroup, record: SAMRecord)
    {
        assert(!readGroup.contains(record))
        if (readGroup.acceptRead(record))
        {
            ++numAcceptedReads
        }
    }

    private fun processCompleteReadGroup(readGroup: ReadGroup)
    {
        if (mIncompleteReadGroups.remove(readGroup.name) != null)
        {
            mIncompleteReadNames.remove(readGroup.name)
            if (readGroup.Reads.size > 2)
            {
                logger.debug("read group size: {}", readGroup.Reads.size)
            }
            writeReadGroup(readGroup)
            ++mNumCompletedGroups
        }
    }

    private fun writeReadGroup(readGroup: ReadGroup)
    {
        for (record in readGroup.allReads)
        {
            addReadDataTableRow(record, readGroup.isComplete())
        }
        for (samRecord in readGroup.allReads)
        {
            mBamFileWriter.addAlignment(samRecord)
        }
    }

    private fun setReadDataTableColumns()
    {
        // add all the columns we need for the CSV
        mReadDataTable.addColumns(
            StringColumn.create("readId"),
            StringColumn.create("chromosome"),
            IntColumn.create("posStart"),
            IntColumn.create("posEnd"),
            StringColumn.create("mateChr"),
            IntColumn.create("matePosStart"),
            BooleanColumn.create("hasTeloContent"),
            StringColumn.create("cigar"),
            IntColumn.create("insertSize"),
            BooleanColumn.create("firstInPair"),
            BooleanColumn.create("unmapped"),
            BooleanColumn.create("mateUnmapped"),
            BooleanColumn.create("isSupplementary"),
            IntColumn.create("flags"),
            StringColumn.create("suppData"),
            BooleanColumn.create("completeFrag"))
    }

    private fun addReadDataTableRow(record: SAMRecord, readGroupIsComplete: Boolean)
    {
        val row = mReadDataTable.appendRow()
        row.setString("readId", record.readName)
        row.setString("chromosome", record.referenceName)
        row.setInt("posStart", record.alignmentStart)
        row.setInt("posEnd", record.alignmentEnd)
        row.setString("mateChr", record.mateReferenceName)
        row.setInt("matePosStart", record.mateAlignmentStart)
        row.setBoolean("hasTeloContent", hasTelomericContent(record.readString))
        row.setString("cigar", record.cigarString)
        row.setInt("insertSize", record.inferredInsertSize)
        row.setBoolean("firstInPair", SamRecordUtils.firstInPair(record))
        row.setBoolean("unmapped", record.readUnmappedFlag)
        row.setBoolean("mateUnmapped", SamRecordUtils.mateUnmapped(record))
        row.setBoolean("isSupplementary", record.supplementaryAlignmentFlag)
        row.setInt("flags", record.flags)
        row.setString("suppData", record.getStringAttribute(SAMTag.SA.name) ?: "")
        row.setBoolean("completeFrag", readGroupIsComplete)
    }

    private fun writeReadDataToTsv()
    {
        if (mReadDataCsvPath == null)
            return

        GZIPOutputStream(java.io.FileOutputStream(mReadDataCsvPath, false)).use { outputStream ->
            try
            {
                val writeOptions = CsvWriteOptions
                    .builder(outputStream)
                    .separator('\t').build()
                mReadDataTable.write().csv(writeOptions)
            } catch (e: java.io.IOException)
            {
                throw IllegalStateException("Could not save to tsv file: " + mReadDataCsvPath + ": " + e.message)
            }
        }
    }

    init
    {
        val samReader = openSamReader(config)
        val telbamPath = config.telbamFile
        mBamFileWriter = SAMFileWriterFactory().makeBAMWriter(samReader.fileHeader, false, java.io.File(telbamPath))
        mReadDataCsvPath = config.tsvFile
        setReadDataTableColumns()
    }
}