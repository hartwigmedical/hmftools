package com.hartwig.hmftools.extensions.csv

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.commons.csv.QuoteMode
import java.io.FileWriter
import kotlin.reflect.KClass
import kotlin.reflect.KParameter
import kotlin.reflect.full.*
import kotlin.reflect.jvm.isAccessible

typealias Prefix = List<String>

object CsvWriter {
    /**
     * Writes a list of CsvData records to file, including header. Will only write fields present in the primary constructor, in their order
     * of declaration and relies on the generated componentN methods of Kotlin data classes to get field values.
     * To avoid naming conflicts between fields with same name defined on different CsvData classes, nested fields are prefixed with their
     * parent field name.
     *
     * e.g:
     *  - the TestCsvData data class will get printed to a csv file with 2 columns, field1 and field2 (field3 will be skipped).
     *  - the NestedCsvData data class will produce a csv file with 3 columns: field0, read.field1 and read.field2.
     * <pre><code>
     *     data class TestCsvData(val field1: String, val field2: String) : CsvData {
     *          val field3: String = ""
     *     }
     *
     *     data class NestedCsvData(val field0: String, val read: TestCsvData): CsvData
     * </code></pre>
     *
     * @param records list of records to print
     * @param location output file location
     * @param format CSV format
     */
    inline fun <reified T : CsvData> write(records: List<T>, location: String, format: CSVFormat) {
        val printer = CSVPrinter(FileWriter(location),
                                 format.withQuoteMode(QuoteMode.MINIMAL).withHeader(*T::class.columns().toTypedArray()))
        printer.printRecords(records.map { T::class.values(it) })
        printer.close()
    }

    inline fun <reified T : CsvData> writeCSV(records: List<T>, location: String, nullString: String = DEFAULT_NULL_STRING) {
        write(records, location, DEFAULT_CSV_FORMAT.withNullString(nullString))
    }

    inline fun <reified T : CsvData> writeTSV(records: List<T>, location: String, nullString: String = DEFAULT_NULL_STRING) {
        write(records, location, DEFAULT_TSV_FORMAT.withNullString(nullString))
    }

    // Generate CsvData header record based on primary constructor parameter list
    fun <T : CsvData> KClass<T>.columns(prefixes: List<String> = listOf()): List<String> {
        if (!this.isData) {
            throw IllegalArgumentException("Cannot write value of type ${this.qualifiedName} to CSV. Must be a data class.")
        }
        return this.primaryConstructor!!.parameters.flatMap { parameterNames(this.qualifiedName, prefixes, it) }
    }

    // Extract CsvData values for primary constructor params
    fun <T : CsvData> KClass<T>.values(csvRecord: T): List<String?> {
        val columnCount = this.columns().size
        return this.declaredFunctions.filter { it.isOperator && it.name.matches(Regex("component[0-9]+")) }
                .sortedBy { it.name.substringAfter("component").toInt() }.flatMap { component ->
                    component.isAccessible = true
                    if (component.returnType.isSubtypeOf(CsvData::class.starProjectedType)) {
                        @Suppress("unchecked_cast")
                        val returnType = component.returnType.classifier as KClass<CsvData>
                        val returnValue = component.call(csvRecord) as CsvData
                        returnType.values(returnValue)
                    } else {
                        listOf(component.call(csvRecord) as String?)
                    }
                }.take(columnCount)
    }

    private fun parameterNames(className: String?, prefix: Prefix, param: KParameter): List<String> {
        return when {
            param.type.isSubtypeOf(String::class.starProjectedType.withNullability(true)) -> listOf(prefixedParamName(prefix, param))
            param.type.isSubtypeOf(CsvData::class.starProjectedType)                      -> {
                @Suppress("unchecked_cast")
                (param.type.classifier as KClass<CsvData>).columns(prefix + param.name!!)
            }
            else                                                                          -> {
                throw IllegalArgumentException("Cannot write value of type $className to CSV. All fields in primary constructor must be " +
                                                       "either String, String? or CsvData")
            }
        }
    }

    private fun prefixedParamName(prefix: Prefix, param: KParameter): String {
        return if (prefix.isEmpty()) {
            param.name!!
        } else {
            prefix.joinToString(".") + "." + param.name!!
        }
    }
}
