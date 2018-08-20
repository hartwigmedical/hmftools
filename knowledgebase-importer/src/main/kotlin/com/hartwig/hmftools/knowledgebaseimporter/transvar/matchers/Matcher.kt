package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

interface Matcher {
    fun matches(string: String): Boolean
    fun contains(string: String): Boolean {
        return string.split("[\\s;/+–]".toRegex()).any { matches(it) }
    }
}
