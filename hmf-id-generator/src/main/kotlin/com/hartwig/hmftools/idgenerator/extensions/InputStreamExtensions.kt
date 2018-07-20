package com.hartwig.hmftools.idgenerator.extensions

import com.hartwig.hmftools.idgenerator.HASH
import com.hartwig.hmftools.idgenerator.HMF_ID
import com.hartwig.hmftools.idgenerator.HmfId
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.io.InputStream
import java.nio.charset.Charset

fun InputStream.readCurrentIds(): Set<HmfId> {
    val parser = CSVParser.parse(this, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader())
    return parser.map { HmfId(it.get(HASH), it.get(HMF_ID).toInt()) }.toSet()
}

