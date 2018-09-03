package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

object TransvarMatcher : Matcher {
    private val nestedMatchers = listOf(TransvarCDnaMatcher, TransvarGDnaMatcher, TransvarProteinMatcher)

    override fun matches(string: String): Boolean {
        return nestedMatchers.any { it.matches(string) }
    }

    override fun contains(string: String, separatorPattern: String) = nestedMatchers.any { it.contains(string) }
}
