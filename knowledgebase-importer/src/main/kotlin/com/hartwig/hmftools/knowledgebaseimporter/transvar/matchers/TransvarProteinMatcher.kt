package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

//MIVO: rules based on:
// - https://github.com/zwdzwd/transvar/blob/v2.4.0.20180701/transvar/mutation.py
object TransvarProteinMatcher : Matcher {
    private val FS_PATTERN = "(p\\.)?(${AminoAcidSymbols.pattern}|[*?])(\\d+)((${AminoAcidSymbols.pattern}|[*?])+)?fs((Ter|[*Xx])(\\d+))?\$"
    private val PATTERN = "(p\\.)?(${AminoAcidSymbols.pattern}|[*?])(\\d+)(_(${AminoAcidSymbols.pattern}|[*?]*)?(\\d+))?(del([^i][A-Za-z*?\\d]*)?)?(ins([A-Za-z*?]+))?(dup)?>?(${AminoAcidSymbols.pattern}|[*?])?(fs((Ter|[*Xx])(\\d+|\\?))?)?(ref([A-Za-zx*]*))?\$"

    override fun matches(string: String): Boolean {
        val input = string.trim()
        return FS_PATTERN.toRegex(RegexOption.IGNORE_CASE).matches(input) || PATTERN.toRegex(RegexOption.IGNORE_CASE).matches(input)
    }
}
