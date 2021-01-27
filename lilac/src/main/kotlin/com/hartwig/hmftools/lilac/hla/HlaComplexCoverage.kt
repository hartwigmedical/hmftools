package com.hartwig.hmftools.lilac.hla

import java.io.File
import kotlin.math.roundToInt

data class HlaComplexCoverage(val uniqueCoverage: Int, val sharedCoverage: Double, val wildCoverage: Double, val alleles: List<HlaAlleleCoverage>) : Comparable<HlaComplexCoverage> {
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

            return HlaComplexCoverage(unique, shared, wild, alleles)
        }

        fun header(): String {
            return "totalCoverage\tuniqueCoverage\tsharedCoverage\twildCoverage\tallele1\tallele2\tallele3\tallele4\tallele5\tallele6"
        }

        fun List<HlaComplexCoverage>.writeToFile(fileName: String) {
            val file = File(fileName)
            file.writeText(header() + "\n" )

            for (coverage in this) {
                file.appendText(coverage.toString() + "\n");
            }
        }
    }


    override fun compareTo(other: HlaComplexCoverage): Int {
        val totalCoverageCompare = totalCoverage.compareTo(other.totalCoverage)
        if (totalCoverageCompare != 0) {
            return totalCoverageCompare
        }

        val uniqueCoverageCompare = uniqueCoverage.compareTo(other.uniqueCoverage)
        if (uniqueCoverageCompare != 0) {
            return uniqueCoverageCompare
        }

        val sharedCoverageCompare = sharedCoverage.compareTo(other.sharedCoverage)
        if (sharedCoverageCompare != 0) {
            return sharedCoverageCompare
        }

        return wildCoverage.compareTo(other.wildCoverage)
    }



    override fun toString(): String {
        return "${totalCoverage.roundToInt()}\t$uniqueCoverage\t${sharedCoverage.roundToInt()}\t${wildCoverage.roundToInt()}\t${alleles.size}\t${alleles.joinToString("\t")}"
    }


}