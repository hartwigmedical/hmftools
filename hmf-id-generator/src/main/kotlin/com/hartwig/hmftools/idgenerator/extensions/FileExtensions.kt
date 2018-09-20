package com.hartwig.hmftools.idgenerator.extensions

import com.hartwig.hmftools.idgenerator.HASH
import com.hartwig.hmftools.idgenerator.HMF_ID
import com.hartwig.hmftools.idgenerator.PATIENT_ID
import com.hartwig.hmftools.idgenerator.anonymizedIds.HashId
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfPatientId
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import java.io.File
import java.nio.charset.Charset

fun File.readOldIds(): Set<HashId> {
    val parser = CSVParser.parse(this, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader())
    return parser.map { HashId(it.get(HASH), it.get(HMF_ID).toInt()) }.toSet()
}

fun File.writeHmfPatientIds(ids: Collection<HmfPatientId>) {
    this.bufferedWriter().use { out ->
        out.write("$HMF_ID,$HASH,$PATIENT_ID")
        out.newLine()
        ids.sortedBy { it.id }.forEach {
            out.write("${it.id},${it.hash},${it.plaintext}")
            out.newLine()
        }
    }
}
