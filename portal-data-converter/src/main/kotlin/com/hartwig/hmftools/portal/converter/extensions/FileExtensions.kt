package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.portal.converter.SampleClinicalData
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.charset.Charset

private val logger = LogManager.getLogger("FileExtensions")

fun File.readSamplesData(): List<SampleClinicalData> {
    val parser = CSVParser.parse(this, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader())
    return parser.map { csvRecord ->
        val sampleData = SampleClinicalData(csvRecord)
        if (sampleData.cpctId == null) {
            logger.warn("Read a CSV record with null CPCT ID")
        }
        sampleData
    }
}
