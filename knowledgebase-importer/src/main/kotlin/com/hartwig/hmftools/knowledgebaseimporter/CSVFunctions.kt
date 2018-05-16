package com.hartwig.hmftools.knowledgebaseimporter

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import java.io.File
import java.io.InputStream
import java.nio.charset.Charset

fun <T> readTSVRecords(fileLocation: String, recordParser: (CSVRecord) -> T): List<T> {
    val format = CSVFormat.TDF.withIgnoreSurroundingSpaces().withFirstRecordAsHeader().withNullString("")
    return readRecords(fileLocation, format, recordParser)
}

fun <T> readTSVRecords(inputStream: InputStream, recordParser: (CSVRecord) -> T): List<T> {
    val format = CSVFormat.TDF.withIgnoreSurroundingSpaces().withFirstRecordAsHeader().withNullString("")
    return readRecords(inputStream, format, recordParser)
}

fun <T> readCSVRecords(fileLocation: String, recordParser: (CSVRecord) -> T): List<T> {
    val format = CSVFormat.DEFAULT.withIgnoreSurroundingSpaces().withFirstRecordAsHeader().withNullString("")
    return readRecords(fileLocation, format, recordParser)
}

fun <T> readCSVRecords(inputStream: InputStream, recordParser: (CSVRecord) -> T): List<T> {
    val format = CSVFormat.DEFAULT.withIgnoreSurroundingSpaces().withFirstRecordAsHeader().withNullString("")
    return readRecords(inputStream, format, recordParser)
}

fun <T> readRecords(fileLocation: String, format: CSVFormat, recordParser: (CSVRecord) -> T): List<T> {
    val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), format)
    return parser.asSequence().map { recordParser(it) }.toList()
}

fun <T> readRecords(inputStream: InputStream, format: CSVFormat, recordParser: (CSVRecord) -> T): List<T> {
    val parser = CSVParser.parse(inputStream, Charset.defaultCharset(), format)
    return parser.asSequence().map { recordParser(it) }.toList()
}
