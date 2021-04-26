package com.hartwig.hmftools.extensions.csv

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import java.io.File
import java.io.InputStream
import java.nio.charset.Charset
import kotlin.reflect.KClass
import kotlin.reflect.full.isSubtypeOf
import kotlin.reflect.full.primaryConstructor
import kotlin.reflect.full.starProjectedType
import kotlin.reflect.full.withNullability

object CsvReader {
    /**
     * Reads CsvData records from a file written with the CsvWriter. Builds CsvData objects by zipping primary constructor params with
     * csv column values (does not take csv column names into account).
     *
     * @param clazz the class of the records to be read. Needs to be a data class.
     * @param fileLocation input file location
     * @param format CSV format
     */
    fun <T : CsvData> read(clazz: KClass<T>, fileLocation: String, format: CSVFormat): List<T> {
        val parser = CSVParser.parse(File(fileLocation), Charset.defaultCharset(), format)
        return parser.asSequence().map { clazz.read(it) }.toList()
    }

    inline fun <reified T : CsvData> readTSV(fileLocation: String, nullString: String? = DEFAULT_NULL_STRING): List<T> =
            read(T::class, fileLocation, DEFAULT_TSV_FORMAT.withNullString(nullString).withFirstRecordAsHeader())

    /**
     * Reads CsvData records from a file. Builds CsvData objects selecting only the CSV columns that match the primary constructor parameters by name
     *
     * @param clazz the class of the records to be read. Needs to be a data class and must not contain nested CsvData fields.
     * @param inputStream input stream
     * @param format CSV format
     */

    fun <T : CsvData> readByName(clazz: KClass<T>, inputStream: InputStream, format: CSVFormat): List<T> {
        val parser = CSVParser.parse(inputStream, Charset.defaultCharset(), format)
        return parser.asSequence().map { clazz.read(it.toMap()) }.toList()
    }

    inline fun <reified T : CsvData> readCSVByName(fileLocation: String, nullString: String? = DEFAULT_NULL_STRING): List<T> =
            readCSVByName(File(fileLocation).inputStream(), nullString)

    inline fun <reified T : CsvData> readCSVByName(inputStream: InputStream, nullString: String? = DEFAULT_NULL_STRING): List<T> =
            readByName(T::class, inputStream, DEFAULT_CSV_FORMAT.withNullString(nullString).withFirstRecordAsHeader())

    private fun <T : CsvData> KClass<T>.read(csvRecord: CSVRecord): T {
        val recordValues = csvRecord.asIterable().toList()
        val (value, remainingValues) = read(recordValues)
        if (remainingValues.isNotEmpty()) {
            throw IllegalArgumentException("Could not correctly map all csv columns on the ${this.qualifiedName} class. Unmapped records: $remainingValues")
        }
        return value
    }

    private fun <T : CsvData> KClass<T>.read(values: List<*>): Pair<T, List<*>> {
        if (!this.isData) {
            throw IllegalArgumentException("Cannot read value of type ${this.qualifiedName} to CSV. Must be a data class.")
        }
        val constructor = this.primaryConstructor!!
        val nestedValues = constructor.parameters.foldIndexed(values) { index, remainingValues, param ->
            when {
                param.type.isSubtypeOf(String::class.starProjectedType.withNullability(true)) -> remainingValues
                param.type.isSubtypeOf(CsvData::class.starProjectedType)                      -> {
                    @Suppress("unchecked_cast")
                    val paramClass = (param.type.classifier as KClass<CsvData>)
                    val fieldsOfParent = remainingValues.take(index)
                    val (nestedValue, unusedValues) = paramClass.read(remainingValues.drop(index))
                    fieldsOfParent + nestedValue + unusedValues
                }
                else                                                                          -> {
                    throw IllegalArgumentException("Cannot read value of type ${this.qualifiedName} from CSV. All fields in primary constructor must be " +
                                                           "either String, String? or CsvData")
                }
            }
        }
        return Pair(this.primaryConstructor!!.call(*nestedValues.take(this.primaryConstructor!!.parameters.size).toTypedArray()),
                    nestedValues.drop(this.primaryConstructor!!.parameters.size))
    }

    private fun <T : CsvData> KClass<T>.read(record: Map<String, String>): T {
        if (!this.isData) {
            throw IllegalArgumentException("Cannot read value of type ${this.qualifiedName} to CSV. Must be a data class.")
        }
        val constructor = this.primaryConstructor!!
        val paramValues = constructor.parameters.map {
            when {
                it.type.isSubtypeOf(String::class.starProjectedType.withNullability(true)) -> {
                    if (!record.containsKey(it.name)) error("Could not find column for field ${it.name} in CSV file.")
                    if (!it.type.isMarkedNullable && record[it.name] == null) error("Found null CSV value for non-nullable field ${it.name}")
                    record[it.name]
                }
                else                                                                       -> {
                    error("Cannot read value of type ${this.qualifiedName} from CSV. All fields in primary constructor must be String or String?")
                }
            }
        }
        return this.primaryConstructor!!.call(*paramValues.toTypedArray())
    }
}
