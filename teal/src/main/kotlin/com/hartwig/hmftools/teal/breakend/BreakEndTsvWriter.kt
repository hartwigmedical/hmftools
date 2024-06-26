package com.hartwig.hmftools.teal.breakend

import org.apache.logging.log4j.LogManager
import tech.tablesaw.api.*
import tech.tablesaw.columns.Column
import tech.tablesaw.io.csv.CsvWriteOptions
import java.io.FileOutputStream
import java.io.IOException
import java.lang.RuntimeException
import java.util.zip.GZIPOutputStream

object BreakEndTsvWriter
{
    private val logger = LogManager.getLogger(javaClass)

    private fun breakEndLogStr(breakEnd: TelomericBreakEndSupport?) : String
    {
        if (breakEnd == null)
        {
            return "missing"
        }

        var logLine = String()

        Fragment.AlignedReadType.values().forEach({ t -> logLine += "${t.code}(${breakEnd.fragmentCount({f -> f.alignedReadType == t})}) " })
        Fragment.PairedReadType.values().forEach({ t -> logLine += "${t.code}(${breakEnd.fragmentCount({f -> f.pairedReadType == t})}) " })

        logLine += "longestSR(${breakEnd.longestSplitReadAlignLength})"
        return logLine
    }

    private fun fragmentTypeToColumn(alignedReadType: Fragment.AlignedReadType,
                                     pairedReadType: Fragment.PairedReadType, prefix: String) : String
    {
        if (alignedReadType == Fragment.AlignedReadType.DISCORDANT_PAIR &&
            pairedReadType == Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC)
        {
            return "${prefix}DPTelNoSR"
        }

        val alignedTypeStr = when (alignedReadType)
        {
            Fragment.AlignedReadType.SPLIT_READ_TELOMERIC -> "SRTel"
            Fragment.AlignedReadType.NOT_SPLIT_READ -> "NotSplit"
            Fragment.AlignedReadType.SPLIT_READ_NOT_TELOMERIC -> "SRNotTel"
            Fragment.AlignedReadType.DISCORDANT_PAIR -> "DP"
        }

        val pairedTypeStr = when (pairedReadType)
        {
            Fragment.PairedReadType.DISCORDANT_PAIR_TELOMERIC -> "DPTel"
            Fragment.PairedReadType.DISCORDANT_PAIR_NOT_TELOMERIC -> "DPNotTel"
            Fragment.PairedReadType.NOT_DISCORDANT -> "NoDP"
        }

        return "${prefix}${alignedTypeStr}${pairedTypeStr}"
    }

    private fun addColumns(breakEndTable: Table)
    {
        val columns = ArrayList<Column<*>>()

        // add all the columns we need for the CSV
        //breakEndTable.addColumns
        //
        columns += listOf(
            StringColumn.create("chromosome"),
            IntColumn.create("position"),
            IntColumn.create("orientation"),
            StringColumn.create("cOrGRich"),
            StringColumn.create("filter"),
            BooleanColumn.create("inTumor"),
            BooleanColumn.create("inGermline"),
            IntColumn.create("distanceToTelomere"),
            IntColumn.create("maxTelomericLength"),
            IntColumn.create("maxAnchorLength"))

        for (prefix in arrayOf("tumor", "germline"))
        {
            // add all the counts columns for the read types
            Fragment.supportFragTypes.forEach({ p -> columns.add(IntColumn.create(fragmentTypeToColumn(p.first, p.second, prefix))) })
            columns.add(IntColumn.create("${prefix}TotalSupport"))
            columns.add(IntColumn.create("${prefix}MAPQ"))
        }

        breakEndTable.addColumns(*columns.toTypedArray())
    }

    private fun populateRow(row: Row, evidence: TelomericBreakEndEvidence)
    {
        val mainBreakEnd = evidence.mainBreakEnd()

        if (mainBreakEnd == null)
        {
            logger.error("invalid breakend evidence, missing both reference and tumor support")
            throw RuntimeException("invalid breakend evidence, missing both reference and tumor support")
        }

        row.setString("chromosome", mainBreakEnd.chromosome)
        row.setInt("position", mainBreakEnd.position)
        row.setInt("orientation", if (mainBreakEnd.type.isRightTelomeric()) 1 else -1)
        if (mainBreakEnd.type == TelomericBreakEndType.RIGHT_G_TELOMERIC ||
            mainBreakEnd.type == TelomericBreakEndType.LEFT_C_TELOMERIC)
        {
            row.setString("cOrGRich", "G")
        }
        else
        {
            row.setString("cOrGRich", "C")
        }
        row.setString("filter", BreakEndEvidenceFilter.getFilterReasons(evidence).joinToString(separator = ","))

        if (mainBreakEnd.distanceToTelomere != null)
            row.setInt("distanceToTelomere", mainBreakEnd.distanceToTelomere!!)

        row.setBoolean("inTumor", evidence.tumorSupport != null)
        row.setBoolean("inGermline", evidence.germlineSupport != null)

        populateData(row, evidence.tumorSupport, "tumor")
        row.setInt("maxTelomericLength", evidence.tumorSupport?.longestTelomereSegment?.length ?: 0)
        row.setInt("maxAnchorLength", evidence.tumorSupport?.longestSplitReadAlignLength ?: 0)

        populateData(row, evidence.germlineSupport, "germline")

        if (logger.isTraceEnabled)
        {
            logger.trace("{} tumor: {}, germline: {}", mainBreakEnd.key, breakEndLogStr(evidence.tumorSupport),
                breakEndLogStr(evidence.germlineSupport))
        }
    }

    fun populateData(row: Row, support: TelomericBreakEndSupport?, prefix: String)
    {
        for ((alignedType, pairedType) in Fragment.supportFragTypes)
        {
            row.setInt(fragmentTypeToColumn(alignedType, pairedType, prefix), support?.fragmentCount(alignedType, pairedType) ?: 0)
        }
        row.setInt("${prefix}TotalSupport", support?.supportFragmentCount() ?: 0)
        row.setInt("${prefix}MAPQ", support?.sumMapQ({ f -> Pair(f.alignedReadType, f.pairedReadType) in Fragment.supportFragTypes }) ?: 0)
    }

    fun writeBreakEnds(csvPath: String, breakEndEvidence: List<TelomericBreakEndEvidence>)
    {
        val breakEndTable = Table.create("TelemereBreakEnds")
        addColumns(breakEndTable)

        for (evidence in breakEndEvidence)
        {
            val row: Row = breakEndTable.appendRow()
            populateRow(row, evidence)
        }

        GZIPOutputStream(FileOutputStream(csvPath, false)).use { outputStream ->
            try
            {
                val writeOptions: CsvWriteOptions = CsvWriteOptions.builder(outputStream).separator('\t').build()
                breakEndTable.write().csv(writeOptions)
            }
            catch (e: IOException)
            {
                throw IllegalStateException("Could not save to tsv file: ${csvPath}: ${e.message}")
            }
        }
    }
}