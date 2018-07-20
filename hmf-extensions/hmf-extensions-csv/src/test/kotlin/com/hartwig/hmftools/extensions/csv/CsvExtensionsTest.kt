package com.hartwig.hmftools.extensions.csv

import io.kotlintest.matchers.shouldBe
import io.kotlintest.matchers.shouldThrow
import io.kotlintest.specs.StringSpec
import java.io.File

class CsvExtensionsTest : StringSpec() {
    data class SimpleDataClass(val value: String) : CsvData
    data class SimpleDataClassNull(val value: String?) : CsvData
    class SimpleClass(val value: String) : CsvData
    data class SimpleListDataClass(val list: List<String>) : CsvData
    data class SimpleNested(val simpleDataClass: SimpleDataClass) : CsvData
    data class SimpleNestedNull(val simpleDataClass: SimpleDataClass?) : CsvData
    data class ComplexNested(val key: String, val simpleNull: SimpleDataClassNull, val nested: SimpleNested) : CsvData

    init {
        "writes and reads simple data class" {
            val createdValue = SimpleDataClass("test")
            readsWhatWasWritten(createdValue)
        }

        "writes and reads simple data class with null value" {
            val createdValue = SimpleDataClassNull(null)
            readsWhatWasWritten(createdValue)
        }

        "does not write any class" {
            val createdValue = SimpleClass("test")
            shouldThrow<IllegalArgumentException> { readsWhatWasWritten(createdValue) }
        }

        "does not write data class with params other than String or CsvData" {
            val createdValue = SimpleListDataClass(listOf())
            shouldThrow<IllegalArgumentException> { readsWhatWasWritten(createdValue) }
        }

        "does not write data class with nested nullable CsvData" {
            val createdValue = SimpleNestedNull(null)
            shouldThrow<IllegalArgumentException> { readsWhatWasWritten(createdValue) }
        }

        "writes and reads nested csvData" {
            val createdValue = SimpleNested(SimpleDataClass("test"))
            readsWhatWasWritten(createdValue)
        }

        "writes and reads complex nested csvData" {
            val createdValue = ComplexNested("key", SimpleDataClassNull(null), SimpleNested(SimpleDataClass("value")))
            readsWhatWasWritten(createdValue)
        }
    }

    private fun createTempFile(): File {
        val tempFile = File.createTempFile("csvExtensions", "testFile")
        tempFile.deleteOnExit()
        return tempFile
    }

    private inline fun <reified T : CsvData> readsWhatWasWritten(value: T) {
        val tempFile = createTempFile()
        CsvWriter.writeTSV(listOf(value), tempFile.path)
        val readValues = CsvReader.readTSV<T>(tempFile.path)
        readValues.size shouldBe 1
        readValues.first() shouldBe value
    }
}
