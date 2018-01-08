package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.portal.converter.SampleClinicalData
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

fun File.readCpctSamples(): Map<String, SampleClinicalData> {
    val parser = CSVParser.parse(this, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader())
    return parser.map { SampleClinicalData(it) }.associateBy { it.sampleId }
}
