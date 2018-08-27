package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.DnaMutations.DNA_DELETION
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.DnaMutations.DNA_DELINS
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.DnaMutations.DNA_DUPLICATION
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.DnaMutations.DNA_INSERTION
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.DnaMutations.DNA_SUBSTITUTION

//MIVO: rules based on:
// https://github.com/zwdzwd/transvar/blob/v2.4.0.20180701/transvar/mutation.py
// https://github.com/biocommons/hgvs/blob/113eeee31b7d4d0a2ca1cfcc0579023992593841/hgvs/_data/hgvs.pymeta
// https://hgvs.readthedocs.io/en/stable/grammar.html
object TransvarCDnaMatcher : Matcher {
    private const val PREFIX = "(c\\.)?"
    private const val NUM = "[\\d]+"
    private const val SIGNED_NUM = "[+-]?$NUM"
    private const val OFFSET = "($SIGNED_NUM)?"
    private const val POSITION = "($SIGNED_NUM$OFFSET|\\*$NUM$OFFSET)"
    private const val INTERVAL = "(${POSITION}_$POSITION|$POSITION)"

    override fun matches(string: String): Boolean {
        val input = string.trim()
        return matchesSubstitution(input) || matchesDeletion(input) || matchesDuplication(input) || matchesInsertion(input) ||
                matchesDelIns(input)
    }

    private fun matchesSubstitution(input: String) = input.matchesPattern("$INTERVAL$DNA_SUBSTITUTION")
    private fun matchesDeletion(input: String) = input.matchesPattern("$INTERVAL$DNA_DELETION")
    private fun matchesDuplication(input: String) = input.matchesPattern("$INTERVAL$DNA_DUPLICATION")
    private fun matchesInsertion(input: String) = input.matchesPattern("$INTERVAL$DNA_INSERTION")
    private fun matchesDelIns(input: String) = input.matchesPattern("$INTERVAL$DNA_DELINS")

    fun type(input: String) = when {
        matchesSubstitution(input) -> SequenceVariantType.SUBSTITUTION
        matchesDeletion(input)     -> SequenceVariantType.DELETION
        matchesDuplication(input)  -> SequenceVariantType.DUPLICATION
        matchesInsertion(input)    -> SequenceVariantType.INSERTION
        matchesDelIns(input)       -> SequenceVariantType.DELINS
        else                       -> SequenceVariantType.OTHER
    }

    private fun String.matchesPattern(pattern: String): Boolean {
        return "$PREFIX$pattern".toRegex(RegexOption.IGNORE_CASE).matches(this)
    }
}
