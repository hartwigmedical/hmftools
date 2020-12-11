package com.hartwig.hmftools.lilac.prot

import java.io.File

object ProteinSequenceFile {

    fun aBoundaries(): List<Int> {
        val boundaries = listOf(24, 114, 206, 298, 337, 348, 365)
        val offsets = listOf(2, 2, 20, 20, 20, 20, 20)

        return boundaries.zip(offsets) { x, y -> x + y }
    }

    fun bBoundaries(): List<Int> {
        val boundaries = listOf(24, 114, 206, 298, 337, 348)
        val offsets = listOf(0, 13, 13, 13, 14, 14)

        return boundaries.zip(offsets) { x, y -> x + y }
    }

    fun cBoundaries(): List<Int> {
        val boundaries = listOf(24, 114, 206, 298, 338, 349, 365)
        val offsets = listOf(0, 5, 24, 24, 30, 36, 36)

        return boundaries.zip(offsets) { x, y -> x + y }
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

    fun readWrappedFile(alignment: String): List<ProteinSequence> {
        val entries = LinkedHashMap<String, ProteinSequence>()
        for (line in File(alignment).readLines()) {
            if (line.startsWith("*", 2)) {
                val splitLine = line.split(" ")
                val allele = splitLine[1].trim()
                val remainder = line.substring(allele.length + 2).replace(" ", "")
                if (entries.containsKey(allele)) {
                    entries[allele] = entries[allele]!!.addProtein(remainder)
                } else {
                    entries[allele] = ProteinSequence(allele, remainder)
                }
            }
        }

        return entries.values.toList()
    }


    fun writeUnwrappedFile(codonBoundaries: List<Int>, input: String, output: String) {

        val originalEntries = readWrappedFile(input)
        val expandedEntries = expand(originalEntries)

        if (originalEntries.isNotEmpty()) {
            val outputFile = File(output)
            outputFile.writeText(header(codonBoundaries))
            for (entry in originalEntries) {
                outputFile.appendText("${entry.contig.padEnd(20, ' ')}\t${entry.proteins}\n")
            }
        }
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
