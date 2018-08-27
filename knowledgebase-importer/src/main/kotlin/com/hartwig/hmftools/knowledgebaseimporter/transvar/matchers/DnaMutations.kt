package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

object DnaMutations {
    private const val DNA = "[GATCN]"
    const val DNA_SUBSTITUTION = "$DNA>$DNA"
    const val DNA_DELETION = "del([\\d]+|$DNA*)"
    const val DNA_DUPLICATION = "dup$DNA*"
    const val DNA_DELINS = "${DNA_DELETION}ins$DNA+"
    const val DNA_INSERTION = "ins$DNA*"
}
