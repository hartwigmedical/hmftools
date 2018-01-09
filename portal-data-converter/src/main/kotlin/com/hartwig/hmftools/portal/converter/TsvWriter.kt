package com.hartwig.hmftools.portal.converter

import com.hartwig.hmftools.portal.converter.records.Record
import com.hartwig.hmftools.portal.converter.records.SampleRecords
import com.hartwig.hmftools.portal.converter.records.donor.Donor
import com.hartwig.hmftools.portal.converter.records.sample.Sample
import com.hartwig.hmftools.portal.converter.records.specimen.Specimen
import com.hartwig.hmftools.portal.converter.records.ssm.SimpleSomaticMutation
import com.hartwig.hmftools.portal.converter.records.ssm.SimpleSomaticMutationMetadata
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.FileWriter
import java.io.IOException
import kotlin.reflect.KClass

object TsvWriter {

    @Throws(IOException::class)
    fun writeSampleRecords(folder: String, records: Collection<SampleRecords>) {
        printRecords(Donor.header, "$folder/donor.", records.map { it.donor }.toSet())
        printRecords(Sample.header, "$folder/sample.", records.map { it.sample })
        printRecords(Specimen.header, "$folder/specimen.", records.map { it.specimen })
    }

    @Throws(IOException::class)
    fun writeSimpleSomaticMutation(folder: String, sampleName: String, simpleSomaticMutations: List<SimpleSomaticMutation>) {
        printRecords(SimpleSomaticMutation.header, "$folder/ssm_p.$sampleName", simpleSomaticMutations)
    }

    @Throws(IOException::class)
    fun writeSomaticMutationMetadata(folder: String, sampleName: String, metadatas: Collection<SimpleSomaticMutationMetadata>) {
        printRecords(SimpleSomaticMutationMetadata.header, "$folder/ssm_m.$sampleName", metadatas)
    }

    private fun <T : Enum<T>> printRecords(headerClass: KClass<T>, fileName: String, records: Collection<Record>) {
        var count = 0
        records.chunked(100000, {
            val printer = createPrinter(headerClass, "$fileName$count.txt")
            printer.printRecords(it.map { it.record() })
            printer.close()
            count++
        })
    }

    @Throws(IOException::class)
    private fun <T : Enum<T>> createPrinter(headerEnum: KClass<T>, fileName: String): CSVPrinter {
        val format = CSVFormat.TDF.withNullString(DEFAULT_VALUE).withHeader(headerEnum.java)
        return CSVPrinter(FileWriter(fileName), format)
    }
}
