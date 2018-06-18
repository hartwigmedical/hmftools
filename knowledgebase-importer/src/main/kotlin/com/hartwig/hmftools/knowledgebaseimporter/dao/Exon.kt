package com.hartwig.hmftools.knowledgebaseimporter.dao

data class Exon(val chromosome: String, val start: Long, val end: Long, val phase: Int, val endPhase: Int) {
    val length = end - start + 1
    val isFirst = phase == -1
    val isLast = endPhase == -1
    val isNonCoding = phase == -1 && endPhase == -1

    // MIVO: offset of the first codon completely contained in this exon:
    // according to http://www.ensembl.org/info/website/glossary.html
    //  - if phase == 1 -> previous intron lands between the 1st and second base of a codon -> skip first 2 bases
    //  - if phase == 2 -> previous intron lands between the second and 3rd base of a codon -> skip first base
    val firstCodonStartOffset = when (phase) {
        1    -> 2
        2    -> 1
        else -> 0
    }
}
