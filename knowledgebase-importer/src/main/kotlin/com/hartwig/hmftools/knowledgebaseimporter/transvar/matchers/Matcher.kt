package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

interface Matcher {
    companion object {
        const val SPLIT_PATTERN = "[\\s;/+–()]"
    }

    fun matches(string: String): Boolean
    fun contains(string: String, separatorPattern: String = SPLIT_PATTERN) = string.split(separatorPattern.toRegex()).any { matches(it) }
}
