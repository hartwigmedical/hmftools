package com.hartwig.hmftools.portal.converter.ids

import java.util.regex.Pattern

fun extractPatientId(sample: String): String? {
    val pattern = Pattern.compile(PATIENT_ID_PATTERN)
    val matcher = pattern.matcher(sample)
    return if (matcher.find()) {
        return matcher.group(0)
    } else {
        null
    }
}
