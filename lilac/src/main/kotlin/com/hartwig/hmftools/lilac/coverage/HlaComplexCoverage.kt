package com.hartwig.hmftools.lilac.coverage

import java.io.File
import kotlin.math.roundToInt

data class HlaComplexCoverage(val uniqueCoverage: Int, val sharedCoverage: Int, val wildCoverage: Int, val alleles: List<HlaAlleleCoverage>) : Comparable<HlaComplexCoverage> {
    val totalCoverage = uniqueCoverage + sharedCoverage + wildCoverage

    companion object {

        fun create(alleles: List<HlaAlleleCoverage>): HlaComplexCoverage {
            var unique = 0
            var shared = 0.0
            var wild = 0.0
            for (coverage in alleles) {
                unique += coverage.uniqueCoverage
                shared += coverage.sharedCoverage
                wild += coverage.wildCoverage
            }

            return HlaComplexCoverage(unique, shared.roundToInt(), wild.roundToInt(), alleles.expand())
        }

        fun header(): String {
            return "totalCoverage\tuniqueCoverage\tsharedCoverage\twildCoverage\ttypes\tallele1\tallele2\tallele3\tallele4\tallele5\tallele6"
        }

        fun List<HlaComplexCoverage>.writeToFile(fileName: String) {
            val file = File(fileName)
            file.writeText(header() + "\n")

            for (coverage in this) {
                file.appendText(coverage.toString() + "\n");
            }
        }

        private fun List<HlaAlleleCoverage>.expand(): List<HlaAlleleCoverage> {
            if (this.size == 6) {
                return this.sortedBy { it.allele }
            }

            val a = this.filter { it.allele.gene == "A" }
            val b = this.filter { it.allele.gene == "B" }
            val c = this.filter { it.allele.gene == "C" }

            return (a.duplicateSingle() + b.duplicateSingle() + c.duplicateSingle()).sortedBy { it.allele }
        }

        private fun List<HlaAlleleCoverage>.duplicateSingle(): List<HlaAlleleCoverage> {
            if (this.size == 1) {
                val single = this[0]
                val duplicate = HlaAlleleCoverage(single.allele, 0, 0.0, 0.0)
                return listOf(single, duplicate)

            }

            return this
        }

    }

    fun homozygousAlleles(): Int {
        val size = alleles.size
        val distinct = alleles.distinct().size
        return size - distinct
    }

    override fun compareTo(other: HlaComplexCoverage): Int {
        val totalCoverageCompare = totalCoverage.compareTo(other.totalCoverage)
        if (totalCoverageCompare != 0) {
            return totalCoverageCompare
        }

        val wildCoverageCompare = wildCoverage.compareTo(other.wildCoverage)
        if (wildCoverageCompare != 0) {
            return -wildCoverageCompare
        }

        val sharedCoverageCompare = sharedCoverage.compareTo(other.sharedCoverage)
        if (sharedCoverageCompare != 0) {
            return sharedCoverageCompare
        }

        return uniqueCoverage.compareTo(other.uniqueCoverage)
    }


    override fun toString(): String {
        val types = alleles.map { it.allele }.toSet()
        return "${totalCoverage}\t$uniqueCoverage\t${sharedCoverage}\t${wildCoverage}\t${types.size}\t${alleles.joinToString("\t")}"
    }


}