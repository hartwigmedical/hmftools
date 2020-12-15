package com.hartwig.hmftools.lilac.prot

import com.hartwig.hmftools.lilac.hla.HlaReferenceSequence
import java.io.File

object UnwrappedProtFile {
    data class ProtEntry(val contig: String, val proteins: String) {
        fun addProtein(more: String): ProtEntry {
            return ProtEntry(contig, proteins + more)
        }

    }

    fun unwrap(fasta: String, alignment: String, output: String) {

        val permittedAlleles = HlaReferenceSequence(File(fasta)).map { x -> x.allele.toString() }.toSet()

        val entries = LinkedHashMap<String, ProtEntry>()
        for (line in File(alignment).readLines()) {
            if (line.startsWith("*", 2)) {
                val splitLine = line.split(" ")
                val allele = splitLine[1].trim()
                if (permittedAlleles.contains(allele)) {
                    val remainder = line.substring(allele.length + 2).replace(" ", "")
                    if (entries.containsKey(allele)) {
                        entries[allele] = entries[allele]!!.addProtein(remainder)
                    } else {
                        entries[allele] = ProtEntry(allele, remainder)
                    }
                }
            }
        }

        val outputFile = File(output)
        outputFile.writeText("")

        for (entry in entries.values) {
            outputFile.appendText("${entry.contig}\t${entry.proteins}\n")

        }

    }
}
