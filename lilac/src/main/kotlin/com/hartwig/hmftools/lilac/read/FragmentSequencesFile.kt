package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.hla.HlaAllele
import java.io.File

object FragmentSequencesFile {

    fun writeFile(filename: String, fragmentSequences: List<FragmentAlleles>) {
        val alleles = (fragmentSequences.flatMap { it.full } union fragmentSequences.flatMap { it.partial }).sortedBy { it.toString() }
        val file = File(filename)
        file.writeText("Fragment\t" + alleles.joinToString("\t") + "\n")

        for (fragment in fragmentSequences) {
            file.appendText(fragment.fragment.id + "\t" + alleles.map { toChar(fragment, it) }.joinToString("\t") + "\n")
        }

    }

    private fun toChar(fragment: FragmentAlleles, allele: HlaAllele): Char {
        if (fragment.full.contains(allele)) {
            return 'T'
        }
        if (fragment.partial.contains(allele)) {
            return 'P'
        }
        return 'F'
    }


}