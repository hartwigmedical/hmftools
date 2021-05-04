package com.hartwig.hmftools.lilackt.qc

import com.hartwig.hmftools.lilackt.SequenceCount
import com.hartwig.hmftools.lilackt.evidence.PhasedEvidence
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci
import org.apache.logging.log4j.LogManager
import kotlin.math.max
import kotlin.math.min

data class HaplotypeQC(val unusedHaplotypes: Int, val unusedHaplotypeMaxSupport: Int, val unusedHaplotypeMaxLength: Int, val unusedHaplotypesPon: Int) {

    companion object {
        private const val MIN_SUPPORT = 7
        private val logger = LogManager.getLogger(this::class.java)
        private val PON_HAPLOTYPES = HaplotypeQC::class.java.getResource("/pon/haplotypes.txt")
                .readText()
                .split("\n")
                .map { Haplotype.fromString(it) }

        private fun Haplotype.inPon(): Boolean {
            return PON_HAPLOTYPES.any { it.contains(this) }
        }

        fun create(minEvidence: Int, winners: Set<HlaSequenceLoci>, evidence: List<PhasedEvidence>, aminoAcidCount: SequenceCount): HaplotypeQC {
            val allUnmatched = evidence
                    .flatMap { it.unmatchedHaplotype(minEvidence, winners, aminoAcidCount) }
                    .sortedByDescending { it.supportingFragments }
                    .distinctBy { x -> (x.startLocus.toString() + x.endLocus.toString() + x.haplotype) }
                    .sortedBy { it.startLocus }

            var pon = 0
            var unusedCount = 0
            var maxSupport = 0
            var maxLength = 0
            for (unmatched in allUnmatched) {

                if (unmatched.supportingFragments >= MIN_SUPPORT) {
                    if (unmatched.inPon()) {
                        pon++
                        logger.info("    UNMATCHED_PON_HAPLTOYPE - $unmatched")
                    } else {
                        maxSupport = max(maxSupport, unmatched.supportingFragments)
                        maxLength = max(maxLength, unmatched.haplotype.length)
                        unusedCount++

                        logger.warn("    UNMATCHED_HAPLTOYPE - $unmatched")
                    }
                }
            }

            return HaplotypeQC(unusedCount, maxSupport, maxLength, pon)
        }

        fun PhasedEvidence.unmatchedHaplotype(minEvidence: Int, winners: Collection<HlaSequenceLoci>, aminoAcidCount: SequenceCount): List<Haplotype> {
            fun consistentWithAny(sequence: String): Boolean {
                return winners.any { it.consistentWith(sequence, *aminoAcidIndices) }
            }

            val unmatched = evidence
                    .filter { !consistentWithAny(it.key) }
                    .filter { it.value >= minEvidence }

            if (unmatched.isEmpty()) {
                return listOf()
            }

            return unmatched.map { Pair(it.key, it.value) }.map { Haplotype.create(this.aminoAcidIndices, it, aminoAcidCount) }
        }
    }

    fun header(): List<String> {
        return listOf("unusedHaplotypes", "unusedHaplotypeMaxSupport", "unusedHaplotypeMaxLength", "unusedHaplotypesPon")
    }

    fun body(): List<String> {
        return listOf(unusedHaplotypes.toString(), unusedHaplotypeMaxSupport.toString(), unusedHaplotypeMaxLength.toString(), unusedHaplotypesPon.toString())
    }

}


data class Haplotype(val startLocus: Int, val endLocus: Int, val supportingFragments: Int, val haplotype: String) {
    companion object {

        fun fromString(line: String): Haplotype {
            val (locus, haplotype) = line.split("\t")
            val locusStart = locus.toInt()
            val locusEnd = locusStart + haplotype.length - 1
            return Haplotype(locus.toInt(), locusEnd, 0, haplotype)
        }

        fun create(aminoAcidIndices: IntArray, evidence: Pair<String, Int>, aminoAcidCount: SequenceCount): Haplotype {
            require(aminoAcidIndices.isNotEmpty())
            val startLoci = aminoAcidIndices.min()!!
            val endLoci = aminoAcidIndices.max()!!
            val sparseHaplotype = evidence.first
            val completeHaplotype = (startLoci..endLoci).joinToString("") { if (aminoAcidIndices.contains(it)) sparseHaplotype[aminoAcidIndices.indexOf(it)].toString() else aminoAcidCount.sequenceAt(it).first() }
            return Haplotype(startLoci, endLoci, evidence.second, completeHaplotype)
        }
    }

    fun contains(unmatched: Haplotype): Boolean {
        if (unmatched.startLocus >= startLocus && unmatched.endLocus <= endLocus) {
            val start = max(unmatched.startLocus, startLocus)
            val end = min(unmatched.endLocus, endLocus)
            for (loci in start..end) {
                if (charAt(loci) != unmatched.charAt(loci)) {
                    return false
                }
            }
            return true
        }

        return false
    }

    private fun charAt(loci: Int): Char {
        require(loci >= startLocus)
        require(loci <= endLocus)

        val index = loci - startLocus
        return haplotype[index]
    }

    override fun toString(): String {
        return "startLocus=$startLocus, endLocus=$endLocus, supportingFragments=$supportingFragments, haplotype=$haplotype"
    }

}