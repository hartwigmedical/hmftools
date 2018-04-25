package com.hartwig.hmftools.knowledgebaseimporter

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import java.io.File
import java.nio.charset.Charset

fun <T> readCSVRecords(fileLocation: String, recordParser: (CSVRecord) -> T): List<T> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), CSVFormat.TDF.withFirstRecordAsHeader().withNullString(""))
    return parser.asSequence().map { recordParser(it) }.toList()
}
