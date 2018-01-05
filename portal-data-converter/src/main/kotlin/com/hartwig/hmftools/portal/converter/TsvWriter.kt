package com.hartwig.hmftools.portal.converter

import com.hartwig.hmftools.portal.converter.donor.Donor
import com.hartwig.hmftools.portal.converter.sample.Sample
import com.hartwig.hmftools.portal.converter.specimen.Specimen
import com.hartwig.hmftools.portal.converter.ssm.SimpleSomaticMutation
import com.hartwig.hmftools.portal.converter.ssm.SimpleSomaticMutationMetadata
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.FileWriter
import java.io.IOException
import kotlin.reflect.KClass

object TsvWriterKt {

    @Throws(IOException::class)
    fun writeDonors(donors: List<Donor>) {
        printRecords(Donor.header, "donor", donors)
    }

    @Throws(IOException::class)
    fun writeSamples(samples: List<Sample>) {
        printRecords(Sample.header, "sample", samples)
    }

    @Throws(IOException::class)
    fun writeSpecimens(specimens: List<Specimen>) {
        printRecords(Specimen.header, "specimen", specimens)
    }

    @Throws(IOException::class)
    fun writeSimpleSomaticMutation(simpleSomaticMutations: List<SimpleSomaticMutation>) {
        printRecords(SimpleSomaticMutation.header, "ssm_p", simpleSomaticMutations)
    }

    @Throws(IOException::class)
    fun writeSomaticMutationMetadata(metadatas: List<SimpleSomaticMutationMetadata>) {
        printRecords(SimpleSomaticMutationMetadata.header, "ssm_m", metadatas)
    }

    private fun <T : Enum<T>> printRecords(headerClass: KClass<T>, fileName: String, records: List<Record>) {
        var count = 0
        records.chunked(100000, {
            val printer = createPrinter(headerClass, "$fileName.$count.txt")
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
