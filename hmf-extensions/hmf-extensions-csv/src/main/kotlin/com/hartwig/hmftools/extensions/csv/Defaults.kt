package com.hartwig.hmftools.extensions.csv

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.QuoteMode

val DEFAULT_CSV_FORMAT: CSVFormat = CSVFormat.DEFAULT.withNullString("!NULL!").withQuoteMode(QuoteMode.MINIMAL).withIgnoreSurroundingSpaces()
val DEFAULT_TSV_FORMAT: CSVFormat = CSVFormat.TDF.withNullString("!NULL!").withQuoteMode(QuoteMode.MINIMAL).withIgnoreSurroundingSpaces()
