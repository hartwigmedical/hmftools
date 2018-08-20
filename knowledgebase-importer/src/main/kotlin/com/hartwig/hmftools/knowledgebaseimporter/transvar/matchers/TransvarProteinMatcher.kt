package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

//MIVO: regex based on: https://github.com/zwdzwd/transvar/blob/v2.4.0.20180701/transvar/mutation.py
object TransvarProteinMatcher : Matcher {
    private val FS_PATTERN = "(p\\.)?(${AminoAcidSymbols.pattern}|[*?])(\\d+)((${AminoAcidSymbols.pattern}|[*?])+)?fs((Ter|[*Xx])(\\d+))?\$"
    private val PATTERN = "(p\\.)?(${AminoAcidSymbols.pattern}|[*?])(\\d+)(_(${AminoAcidSymbols.pattern}|[*?]*)?(\\d+))?(del([^i][A-Za-z*?\\d]*)?)?(ins([A-Za-z*?]+))?>?([A-Za-z*?]+)?(fs((Ter|[*Xx])(\\d+))?)?(ref([A-Za-zx*]*))?\$"

    override fun matches(string: String): Boolean {
        val trimmedInput = string.trim()
        return FS_PATTERN.toRegex().matches(trimmedInput) || PATTERN.toRegex().matches(trimmedInput)
    }
}
