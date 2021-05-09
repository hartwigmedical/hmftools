package com.hartwig.hmftools.lilackt.cna

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage
import com.hartwig.hmftools.lilackt.hla.HlaAllele

data class HlaCopyNumber(val allele: HlaAllele, val copyNumber: Double) {

    companion object {
        private val genes = listOf("HLA-A", "HLA-B", "HLA-C")

        fun alleleCopyNumber(winners: List<HlaAllele>): List<HlaCopyNumber> {
            return winners.map { HlaCopyNumber(it, 0.0) }
        }

        fun alleleCopyNumber(winners: List<HlaAllele>, geneCopyNumberFile: String, tumorCoverage: HlaComplexCoverage): List<HlaCopyNumber> {
            if (geneCopyNumberFile.isEmpty() || tumorCoverage.alleleCoverage.isEmpty()) {
                return winners.map { HlaCopyNumber(it, 0.0) }
            }

            val coverage = tumorCoverage.alleleCoverage
            require(coverage.size == 6)
            val geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberFile).filter { it.gene() in genes}.sortedBy { it.gene() }
            val aCopyNumber = alleleCopyNumber(geneCopyNumbers[0], coverage.filter { it.allele.gene == "A" })
            val bCopyNumber = alleleCopyNumber(geneCopyNumbers[1], coverage.filter { it.allele.gene == "B" })
            val cCopyNumber = alleleCopyNumber(geneCopyNumbers[2], coverage.filter { it.allele.gene == "C" })

            return (aCopyNumber + bCopyNumber + cCopyNumber).sortedBy { it.allele }
        }

        private fun alleleCopyNumber(geneCopyNumber: GeneCopyNumber, alleleCoverage: List<HlaAlleleCoverage>): List<HlaCopyNumber> {
            require(alleleCoverage.size == 2)

            val minor = geneCopyNumber.minMinorAlleleCopyNumber()
            val major = geneCopyNumber.minCopyNumber() - minor

            return if (alleleCoverage[0].totalCoverage >= alleleCoverage[1].totalCoverage) {
               listOf(HlaCopyNumber(alleleCoverage[0].allele, major), HlaCopyNumber(alleleCoverage[1].allele, minor))
            } else {
                listOf(HlaCopyNumber(alleleCoverage[0].allele, minor), HlaCopyNumber(alleleCoverage[1].allele, major))
            }
        }

    }

}
