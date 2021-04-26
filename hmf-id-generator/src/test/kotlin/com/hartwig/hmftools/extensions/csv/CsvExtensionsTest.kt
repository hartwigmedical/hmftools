package com.hartwig.hmftools.extensions.csv

import com.google.common.io.Resources
import io.kotlintest.matchers.shouldBe
import io.kotlintest.matchers.shouldThrow
import io.kotlintest.specs.StringSpec
import java.io.File

class CsvExtensionsTest : StringSpec() {

    private val fileLocation = Resources.getResource("csv/test.csv").path

    data class SimpleDataClass(val value: String) : CsvData
    data class SimpleDataClassNull(val value: String?) : CsvData
    class SimpleClass(val value: String) : CsvData
    data class SimpleListDataClass(val list: List<String>) : CsvData
    data class SimpleNested(val simpleDataClass: SimpleDataClass) : CsvData
    data class SimpleNestedNull(val simpleDataClass: SimpleDataClass?) : CsvData
    data class ComplexNested(val key: String, val simpleNull: SimpleDataClassNull, val nested: SimpleNested) : CsvData
    data class Col1(val col1: String) : CsvData
    data class Col12(val col1: String, val col2: String) : CsvData
    data class Col3(val col3: String) : CsvData
    data class ColNull(val colNull: String?) : CsvData
    data class ColNonNull(val colNull: String) : CsvData
    data class Col4(val col4: String) : CsvData

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

        "reads by name column 1 from file" {
            val records = CsvReader.readCSVByName<Col1>(fileLocation)
            records[0] shouldBe Col1("v1")
        }

        "reads by name column 1 and 2 from file" {
            val records = CsvReader.readCSVByName<Col12>(fileLocation)
            records[0] shouldBe Col12("v1", "v2")
        }

        "reads by name only column 3 from file" {
            val records = CsvReader.readCSVByName<Col3>(fileLocation)
            records[0] shouldBe Col3("v3")
        }

        "reads by name column with null value from file" {
            val records = CsvReader.readCSVByName<ColNull>(fileLocation)
            records[0] shouldBe ColNull(null)
        }

        "throws on column with null value for non-nullable class field" {
            shouldThrow<IllegalStateException> { CsvReader.readCSVByName<ColNonNull>(fileLocation) }
        }

        "throws on column class field that does not match any column name" {
            shouldThrow<IllegalStateException> { CsvReader.readCSVByName<Col4>(fileLocation) }
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
