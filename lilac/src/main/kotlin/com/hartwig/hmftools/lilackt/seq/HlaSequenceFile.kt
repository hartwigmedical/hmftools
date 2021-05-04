package com.hartwig.hmftools.lilackt.seq

import com.hartwig.hmftools.lilackt.hla.HlaAllele
import java.io.File

object HlaSequenceFile {

    fun readFile(filename: String): List<HlaSequence> {
        val entries = LinkedHashMap<String, HlaSequence>()
        for (line in File(filename).readLines()) {
            val trimmed = line.trim()

            if (trimmed.startsWith("*", 1)) {
                val split = trimmed.split(" ")
                val allele = split[0].trim()
                val alleleIndex = line.indexOf(allele)
                val remainder = line.substring(alleleIndex + allele.length).trim().replace(" ", "")
                if (entries.containsKey(allele)) {
                    entries[allele] = entries[allele]!!.copyWithAdditionalSequence(remainder)
                } else {
                    entries[allele] = HlaSequence(HlaAllele(allele), remainder)
                }
            }
        }

        return entries.values.toList()
    }

    fun List<HlaSequence>.reduceToFourDigit(): List<HlaSequence> {
        return reduce { it.asFourDigit() }
    }

    fun List<HlaSequence>.reduceToSixDigit(): List<HlaSequence> {
        return reduce { it.asSixDigit() }
    }

    private fun List<HlaSequence>.reduce(transform: (HlaAllele) -> HlaAllele): List<HlaSequence> {
        val resultMap = LinkedHashMap<HlaAllele, HlaSequence>()

        for (sequence in this) {
            val reducedAllele = transform(sequence.allele)
            if (!resultMap.containsKey(reducedAllele)) {
                resultMap[reducedAllele] = HlaSequence(reducedAllele, sequence.rawSequence)
            }
        }

        return resultMap.values.toList()
    }
}