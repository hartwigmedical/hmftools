package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.portal.converter.SampleClinicalData
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

fun File.readSamplesData(): List<SampleClinicalData> {
    val parser = CSVParser.parse(this, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader())
    return parser.map { SampleClinicalData(it) }
}
