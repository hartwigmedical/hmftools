package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

interface Matcher {
    companion object {
        const val SPLIT_PATTERN = "[\\s;/+–()]"
    }

    fun matches(string: String): Boolean
    fun contains(string: String): Boolean {
        return string.split(SPLIT_PATTERN.toRegex()).any { matches(it) }
    }
}
