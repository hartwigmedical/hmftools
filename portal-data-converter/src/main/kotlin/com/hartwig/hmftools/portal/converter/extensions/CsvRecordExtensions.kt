package com.hartwig.hmftools.portal.converter.extensions

import org.apache.commons.csv.CSVRecord

fun CSVRecord.getNullable(field: String): String? {
    val fieldValue = this.get(field)
    return if (fieldValue.isNullOrBlank() || fieldValue.toLowerCase() == "null") {
        null
    } else {
        fieldValue
    }
}
