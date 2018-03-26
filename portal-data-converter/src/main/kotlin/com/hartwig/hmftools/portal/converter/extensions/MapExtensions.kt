package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.portal.converter.DEFAULT_VALUE

inline fun <reified T : Enum<T>> Map<T, String>.toRecord(): List<String> {
    return enumValues<T>().map { this[it] ?: DEFAULT_VALUE }
}
