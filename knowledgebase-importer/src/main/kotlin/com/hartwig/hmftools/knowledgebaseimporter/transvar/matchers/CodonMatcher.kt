package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

//MIVO: Source: https://github.com/zwdzwd/transvar/blob/v2.4.0.20180701/transvar/mutation.py
object CodonMatcher : Matcher {
    private val AA = "(${AminoAcidSymbols.pattern})"
    private val PATTERN = "$AA(\\d+)\$"

    override fun matches(string: String): Boolean {
        val input = string.trim()
        return PATTERN.toRegex(RegexOption.IGNORE_CASE).matches(input)
    }
}
