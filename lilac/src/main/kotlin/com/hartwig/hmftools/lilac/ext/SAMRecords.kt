package com.hartwig.hmftools.lilac.ext

import com.hartwig.hmftools.common.codon.Codons
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.SAMRecord

fun SAMRecord.containsIndel(): Boolean {
    for (cigarElement in this.cigar) {
        if (cigarElement.operator == CigarOperator.I || cigarElement.operator == CigarOperator.D) {
            return true
        }
    }
    return false
}


fun SAMRecord.toAminoAcids(): List<String> {
    val result = mutableListOf<String>()
    val forward = this.readString
    val reverse = forward.reversed()

    result.add(Codons.aminoAcids(forward))
    result.add(Codons.aminoAcids(forward.substring(1)))
    result.add(Codons.aminoAcids(forward.substring(2)))

    result.add(Codons.aminoAcids(reverse))
    result.add(Codons.aminoAcids(reverse.substring(1)))
    result.add(Codons.aminoAcids(reverse.substring(2)))

    return result
}