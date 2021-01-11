package com.hartwig.hmftools.lilac.seq

import org.apache.logging.log4j.LogManager
import java.io.File

object HlaSequenceFile {
    val logger = LogManager.getLogger(this::class.java)

    fun readFile(alignment: String): List<HlaSequence> {
        val entries = LinkedHashMap<String, HlaSequence>()
        for (line in File(alignment).readLines()) {
            val trimmed = line.trim()

            if (trimmed.startsWith("*", 1)) {
                val split = trimmed.split(" ")
                val allele = split[0].trim()
                val alleleIndex = line.indexOf(allele)
                val remainder = line.substring(alleleIndex + allele.length).trim().replace(" ", "")
                if (entries.containsKey(allele)) {
                    entries[allele] = entries[allele]!!.copyWithAdditionalSequence(remainder)
                } else {
                    entries[allele] = HlaSequence(allele, remainder)
                }
            }
        }

        return entries.values.toList()
    }

    fun List<HlaSequence>.inflate(): List<HlaSequence> {
        val template = this[0].rawSequence
        return this.map {it.inflate(template) }
    }

    fun List<HlaSequence>.deflate(): List<HlaSequence> {
        if (this.isEmpty()) {
            return this
        }

        val template = this[0].sequence
        return listOf(this[0]) + this.drop(1).map { it.deflate(template) }
    }

    fun List<HlaSequence>.reduceToFirstFourDigits(): List<HlaSequence> {
        val resultMap = LinkedHashMap<String, HlaSequence>()

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

    fun wipeFile(outputFileName: String) {
        val outputFile = File(outputFileName)
        outputFile.writeText("")
    }

    fun writeBoundary(boundaries: Collection<Int>, outputFileName: String) {
        val outputFile = File(outputFileName)
        outputFile.appendText("Boundary".padEnd(20, ' ') + "\t")
        for (i in 0.. boundaries.max()!!) {
            if (i in boundaries) {
                outputFile.appendText("|")
            } else {
                outputFile.appendText(" ")
            }
        }

        outputFile.appendText("\n")
    }

    fun appendFile(outputFileName: String, sequences: List<HlaSequence>) {
        if (sequences.isNotEmpty()) {
            val outputFile = File(outputFileName)
            for (sequence in sequences) {
                outputFile.appendText(sequence.toString() + "\n")
            }
        }
    }

    fun writeFile(outputFileName: String, sequences: List<HlaSequence>) {
        if (sequences.isNotEmpty()) {
            val outputFile = File(outputFileName)
            outputFile.writeText("")

            for (sequence in sequences) {
                outputFile.appendText(sequence.toString() + "\n")
            }
        }
    }

}