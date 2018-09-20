package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.HgvsVariantType

//MIVO: rules based on:
// https://github.com/zwdzwd/transvar/blob/v2.4.0.20180701/transvar/mutation.py
// https://github.com/biocommons/hgvs/blob/113eeee31b7d4d0a2ca1cfcc0579023992593841/hgvs/_data/hgvs.pymeta
// https://hgvs.readthedocs.io/en/stable/grammar.html
object TransvarGDnaMatcher : Matcher {
    private const val PREFIX = "(g\\.)?"
    private const val NUM = "[\\d]+"
    private const val POSITION = "($NUM)"
    private const val INTERVAL = "(${POSITION}_$POSITION|$POSITION)"

    override fun matches(string: String): Boolean {
        val input = string.trim()
        return matchesSubstitution(input) || matchesDeletion(input) || matchesDuplication(input) || matchesInsertion(input) ||
                matchesDelIns(input)
    }

    private fun matchesSubstitution(input: String) = input.matchesPattern("$INTERVAL${DnaMutations.DNA_SUBSTITUTION}")
    private fun matchesDeletion(input: String) = input.matchesPattern("$INTERVAL${DnaMutations.DNA_DELETION}")
    private fun matchesDuplication(input: String) = input.matchesPattern("$INTERVAL${DnaMutations.DNA_DUPLICATION}")
    private fun matchesInsertion(input: String) = input.matchesPattern("$INTERVAL${DnaMutations.DNA_INSERTION}")
    private fun matchesDelIns(input: String) = input.matchesPattern("$INTERVAL${DnaMutations.DNA_DELINS}")

    private fun String.matchesPattern(pattern: String): Boolean {
        return "$PREFIX$pattern".toRegex(RegexOption.IGNORE_CASE).matches(this)
    }

    fun type(input: String) = when {
        matchesSubstitution(input) -> HgvsVariantType.SUBSTITUTION
        matchesDeletion(input)     -> HgvsVariantType.DELETION
        matchesDuplication(input)  -> HgvsVariantType.DUPLICATION
        matchesInsertion(input)    -> HgvsVariantType.INSERTION
        matchesDelIns(input)       -> HgvsVariantType.DELINS
        else                       -> HgvsVariantType.OTHER
    }
}
