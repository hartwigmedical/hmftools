package com.hartwig.hmftools.extensions.csv

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode

const val DEFAULT_NULL_STRING = "!NULL!"
val DEFAULT_CSV_FORMAT: CSVFormat = CSVFormat.DEFAULT.withQuoteMode(QuoteMode.MINIMAL).withIgnoreSurroundingSpaces()
val DEFAULT_TSV_FORMAT: CSVFormat = CSVFormat.TDF.withQuoteMode(QuoteMode.MINIMAL).withIgnoreSurroundingSpaces()
