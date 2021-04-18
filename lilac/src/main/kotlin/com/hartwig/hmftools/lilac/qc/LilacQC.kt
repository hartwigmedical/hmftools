package com.hartwig.hmftools.lilac.qc

import java.io.File


data class LilacQC(val status: Set<LilacQCStatus>, val aminoAcidQC: AminoAcidQC, val bamQC: BamQC, val coverageQC: CoverageQC, val haplotypeQC: HaplotypeQC, val somaticVariantQC: SomaticVariantQC) {

    companion object {
        fun create(aminoAcidQC: AminoAcidQC, bamQC: BamQC, coverageQC: CoverageQC, haplotypeQC: HaplotypeQC, somaticVariantQC: SomaticVariantQC): LilacQC {
            val status = mutableSetOf<LilacQCStatus>()

            if (haplotypeQC.unusedHaplotypes > 0) {
                status.add(LilacQCStatus.WARN_UNMATCHED_HAPLOTYPE)
            }

            if (coverageQC.aTypes == 0 || coverageQC.bTypes == 0 || coverageQC.cTypes == 0 ) {
                status.add(LilacQCStatus.WARN_UNMATCHED_TYPE)
            }

            if (bamQC.discardedIndelFragments > 0) {
                status.add(LilacQCStatus.WARN_UNMATCHED_INDEL)
            }

            if (somaticVariantQC.unmatchedVariants()) {
                status.add(LilacQCStatus.WARN_UNMATCHED_SOMATIC_VARIANT)
            }

            if (coverageQC.percentWildcard > 0) {
                status.add(LilacQCStatus.WARN_WILDCARD_MATCH)
            }

            if (status.isEmpty()) {
                status.add(LilacQCStatus.PASS)
            }

            return LilacQC(status, aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC)
        }
    }

    fun header(): List<String> {
        return listOf("status") + aminoAcidQC.header() + bamQC.header() + coverageQC.header() + haplotypeQC.header() + somaticVariantQC.header()
    }

    fun body(): List<String> {
        return listOf(status.joinToString(",")) + aminoAcidQC.body() + bamQC.body() + coverageQC.body() + haplotypeQC.body() + somaticVariantQC.body()
    }

    fun writefile(fileName: String) {
        val header = header().joinToString(separator = "\t")
        val body = body().joinToString(separator = "\t")

        val file = File(fileName)
        file.writeText(header + "\n")
        file.appendText(body + "\n")
    }
}

