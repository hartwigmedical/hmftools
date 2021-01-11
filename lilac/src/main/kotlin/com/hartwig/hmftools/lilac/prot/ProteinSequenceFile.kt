package com.hartwig.hmftools.lilac.prot

import com.hartwig.hmftools.lilac.hla.HlaAllele
import org.apache.logging.log4j.LogManager
import java.io.File

object ProteinSequenceFile {

    val logger = LogManager.getLogger(this::class.java)

    fun readWrappedFile(alignment: String): List<ProteinSequence> {
        return readFile(alignment).inflate().firstFourDigits()
    }

    fun writeUnwrappedFile(codonBoundaries: List<Int>, input: String, output: String) {

        val indistinguisable = setOf(
                HlaAllele("A*01:01:01:01"),
                HlaAllele("A*01:335"),
                HlaAllele("A*01:338"),
                HlaAllele("A*02:01:01:01"),
                HlaAllele("A*02:96"),
                HlaAllele("A*02:498"),
                HlaAllele("A*02:716"),
                HlaAllele("A*02:844"),
                HlaAllele("A*02:870"),
                HlaAllele("A*02:891")

        )

//        val specificSequences = setOf(HlaAllele("C*01:02:01:01"), HlaAllele("C*07:01:01:01"), HlaAllele("C*07:877"), HlaAllele("C*07:879"), HlaAllele("C*07:882"))
        val bComplex = setOf(HlaAllele("B*07:02:01:01"), HlaAllele("B*08:01:01:01"), HlaAllele("B*56:01:01:01"), HlaAllele("B*56:68"))
        val bComplexLong = setOf(
                HlaAllele("B*07:02:01:01"),
                HlaAllele("B*08:01:01:01"),
                HlaAllele("B*08:20:01"),
                HlaAllele("B*08:26:03"),
                HlaAllele("B*08:39"),
                HlaAllele("B*08:132"),
                HlaAllele("B*08:134"),
                HlaAllele("B*08:151"),
                HlaAllele("B*08:221"),
                HlaAllele("B*08:225"),
                HlaAllele("B*08:230"),
                HlaAllele("B*08:233"),
                HlaAllele("B*08:248"),
                HlaAllele("B*42:21"),
                HlaAllele("B*55:52"),
                HlaAllele("B*56:01:01:01"),
                HlaAllele("B*56:68"))


        val originalEntries = readFile(input)
        val fourDigitEntries = originalEntries.firstFourDigits()
//                .inflate()
                .filter { it.allele in indistinguisable }

        if (fourDigitEntries.isNotEmpty()) {
            val outputFile = File(output)
            outputFile.writeText(header(codonBoundaries))
            for (entry in fourDigitEntries) {
                val protein = entry.protein
                outputFile.appendText("${entry.contig.padEnd(20, ' ')}\t${protein}\n")
            }
        }
    }

    fun readFile(alignment: String): List<ProteinSequence> {
        val entries = LinkedHashMap<String, ProteinSequence>()
        for (line in File(alignment).readLines()) {
            if (line.startsWith("*", 2)) {
                val splitLine = line.split(" ")
                val allele = splitLine[1].trim()
                val remainder = line.substring(allele.length + 2).replace(" ", "")
                if (entries.containsKey(allele)) {
                    entries[allele] = entries[allele]!!.copyWithAdditionalProtein(remainder)
                } else {
                    entries[allele] = ProteinSequence(allele, remainder)
                }
            }
        }

        return entries.values.toList()
    }


    private fun header(boundaryIndices: List<Int>): String {

        val builder = StringBuilder()
                .append("ExonBoundaries".padEnd(20, ' '))
                .append("\t")

        var previousIndex = -1
        for (i in boundaryIndices.indices) {
            val currentIndex = boundaryIndices[i]
            val spaces = currentIndex - previousIndex - 1
            if (spaces > 0) {
                builder.append(" ".repeat(spaces))
            }
            builder.append("|")
            previousIndex = currentIndex
        }


        return builder.append("\n").toString()
    }


    fun List<ProteinSequence>.firstFourDigits(): List<ProteinSequence> {
        val resultMap = LinkedHashMap<String, ProteinSequence>()

        for (sequence in this) {
            val fourDigitName = sequence.allele.fourDigitName()
            if (!resultMap.containsKey(fourDigitName)) {
                resultMap[fourDigitName] = sequence
            } else {
                val existing = resultMap[fourDigitName]!!
                if (existing.length < sequence.length) {
                    logger.debug("Replacing ${existing.allele} with ${sequence.allele}")
                    resultMap[fourDigitName] = sequence
                }
            }

        }

        return resultMap.values.toList()
    }

    fun List<ProteinSequence>.inflate(): List<ProteinSequence> {
        val template = this[0].protein
        return this.map { ProteinSequence(it.contig, inflate(template, it.protein)) }
    }

    private fun inflate(template: String, sequence: String): String {
        val joiner = StringBuilder()
        for (i in sequence.indices) {
            if (sequence[i] == '-') {
                joiner.append(template[i])
            } else {
                joiner.append(sequence[i])
            }
        }

        return joiner.toString()
    }

}
