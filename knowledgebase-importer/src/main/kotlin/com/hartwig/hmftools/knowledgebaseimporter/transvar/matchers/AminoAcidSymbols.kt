package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

object AminoAcidSymbols {
    //MIVO: based on https://github.com/zwdzwd/transvar/blob/41df7ceab8e0ad1881d703674a54b2d540cb3c40/transvar/utils.py#L101

    private val oneToThree = mapOf("G" to "Gly",
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
                                   "T" to "Thr",
                                   "Z" to "Glx",
                                   "B" to "Asx")

    private val threeToOne = oneToThree.map { (letter, letters) -> Pair(letters, letter) }.toMap()
    private val oneLetterPattern = "[${oneToThree.map { it.key }.joinToString("")}]"
    private val threeLetterPattern = "(${threeToOne.map { it.key }.joinToString("|")})"
    val pattern = "($oneLetterPattern|$threeLetterPattern)"
}
