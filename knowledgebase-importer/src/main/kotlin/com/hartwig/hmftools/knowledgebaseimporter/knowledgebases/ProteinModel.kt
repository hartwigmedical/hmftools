package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

object ProteinModel {
    private val proteinCodes = mapOf("G" to "Gly",
                                     "A" to "Ala",
                                     "L" to "Leu",
                                     "M" to "Met",
                                     "F" to "Phe",
                                     "W" to "Trp",
                                     "K" to "Lys",
                                     "Q" to "Gln",
                                     "E" to "Glu",
                                     "S" to "Ser",
                                     "P" to "Pro",
                                     "V" to "Val",
                                     "I" to "Ile",
                                     "C" to "Cys",
                                     "Y" to "Tyr",
                                     "H" to "His",
                                     "R" to "Arg",
                                     "N" to "Asn",
                                     "D" to "Asp",
                                     "T" to "Thr")

    private val singleLetterPattern = "[${proteinCodes.map { it.key }.joinToString()}]"
    private val triLetterPattern = "(?:${proteinCodes.map { it.value }.joinToString("|")})"
    private const val codonNumberPattern = "[1-9][0-9]*"
    private val singleLetterMutationPattern = Regex("$singleLetterPattern$codonNumberPattern$singleLetterPattern", RegexOption.IGNORE_CASE)
    private val triLetterMutationPattern = Regex("$triLetterPattern$codonNumberPattern$triLetterPattern", RegexOption.IGNORE_CASE)
    private val singleLetterCodon = Regex("$singleLetterPattern$codonNumberPattern", RegexOption.IGNORE_CASE)
    private val triLetterCodon = Regex("$triLetterPattern$codonNumberPattern", RegexOption.IGNORE_CASE)

    fun isProteinAlteration(mutation: String): Boolean {
        return singleLetterMutationPattern.matches(mutation) || triLetterMutationPattern.matches(mutation)
    }

    fun isCodonAlteration(mutation: String): Boolean {
        return singleLetterCodon.matches(mutation) || triLetterCodon.matches(mutation)
    }

    fun hasProteinAlteration(mutation: String): Boolean {
        return singleLetterMutationPattern.containsMatchIn(mutation) || triLetterMutationPattern.containsMatchIn(mutation)
    }

    fun hasCodonAlteration(mutation: String, strict: Boolean = true): Boolean {
        return (singleLetterCodon.containsMatchIn(mutation) || triLetterCodon.containsMatchIn(mutation)) &&
                !(strict && hasProteinAlteration(mutation))
    }
}
