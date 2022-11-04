package com.hartwig.hmftools.cider.primer

import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedReader
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import java.io.IOException

object PrimerTsvFile
{
    private enum class Column
    {
        gene,
        probeName,
        sequence,
        vj
    }

    private val sLogger = LogManager.getLogger(VDJSequenceTsvWriter::class.java)

    @Throws(IOException::class)
    fun load(path: String): List<Primer>
    {
        val primerList: MutableList<Primer> = ArrayList()

        createBufferedReader(path).use { reader ->
            val format = CSVFormat.Builder.create()
                .setDelimiter('\t')
                .setRecordSeparator('\n')
                .setHeader().setSkipHeaderRecord(true) // use first line header as column names
                .build()
            val records: Iterable<CSVRecord> = format.parse(reader)
            for (record in records)
            {
                val geneName = record[Column.gene]
                val probeName = record[Column.probeName]
                val sequence = record[Column.sequence]
                val vj = VJ.valueOf(record[Column.vj])
                primerList.add(Primer(geneName, probeName, sequence, vj))
            }
        }
        return primerList
    }
}