package com.hartwig.hmftools.extensions.csv

/**
 * Represents a CSV record. This interface should be extended only by data classes that have only String/String? or nested CsvData fields
 * in their primary constructor. Such data classes can be easily written and read from files using the CsvWriter and CsvReader objects.
 */
interface CsvData
