package com.hartwig.hmftools.lilac.prot

import org.apache.logging.log4j.LogManager
import java.io.File

object ProteinSequenceFile {

    val logger = LogManager.getLogger(this::class.java)


    fun readWrappedFile(alignment: String): List<ProteinSequence> {
        return reduceToFourDigit(expand(readFile(alignment)))
    }

    fun writeUnwrappedFile(codonBoundaries: List<Int>, input: String, output: String) {
        val originalEntries = readFile(input)
        val fourDigitEntries = reduceToFourDigit(originalEntries)

        if (fourDigitEntries.isNotEmpty()) {
            val outputFile = File(output)
            outputFile.writeText(header(codonBoundaries))
            for (entry in fourDigitEntries) {
                outputFile.appendText("${entry.contig.padEnd(20, ' ')}\t${entry.proteins}\n")
            }
        }
    }

    private fun readFile(alignment: String): List<ProteinSequence> {
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


    private fun reduceToFourDigit(sequences: List<ProteinSequence>): List<ProteinSequence> {
        val map = LinkedHashMap<String, ProteinSequence>()

        for (sequence in sequences) {
            val fourDigitName = sequence.allele.fourDigitName()
            if (!map.containsKey(fourDigitName)) {
                map[fourDigitName] = sequence
            } else {
                val existing = map[fourDigitName]!!
                if (existing.length < sequence.length) {
                    logger.info("Replacing ${existing.allele} with ${sequence.allele}")
                    map[fourDigitName] = sequence
                }
            }

        }

        return map.values.toList()
    }

    private fun expand(sequences: List<ProteinSequence>): List<ProteinSequence> {
        val template = sequences[0].proteins
        return sequences.map { ProteinSequence(it.contig, expand(template, it.proteins)) }
    }

    private fun expand(template: String, sequence: String): String {
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
