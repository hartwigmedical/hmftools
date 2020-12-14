package com.hartwig.hmftools.lilac.prot

import org.apache.logging.log4j.LogManager
import java.io.File

object ProteinSequenceFile {

    val logger = LogManager.getLogger(this::class.java)

    fun readWrappedFile(alignment: String): List<ProteinSequence> {
        return readFile(alignment).inflate().firstFourDigits()
    }

    fun writeUnwrappedFile(codonBoundaries: List<Int>, input: String, output: String) {
        val originalEntries = readFile(input)
        val fourDigitEntries = originalEntries.firstFourDigits()//.filter { it.allele.fourDigitName() == "A*01:01" || it.allele.fourDigitName() == "A*11:01" }

        if (fourDigitEntries.isNotEmpty()) {
            val outputFile = File(output)
            outputFile.writeText(header(codonBoundaries))
            for (entry in fourDigitEntries) {
                outputFile.appendText("${entry.contig.padEnd(20, ' ')}\t${entry.proteins}\n")
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
        val map = LinkedHashMap<String, ProteinSequence>()

        for (sequence in this) {
            val fourDigitName = sequence.allele.fourDigitName()
            if (!map.containsKey(fourDigitName)) {
                map[fourDigitName] = sequence
            } else {
                val existing = map[fourDigitName]!!
                if (existing.length < sequence.length) {
                    logger.debug("Replacing ${existing.allele} with ${sequence.allele}")
                    map[fourDigitName] = sequence
                }
            }

        }

        return map.values.toList()
    }

    fun List<ProteinSequence>.inflate(): List<ProteinSequence> {
        val template = this[0].proteins
        return this.map { ProteinSequence(it.contig, inflate(template, it.proteins)) }
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
