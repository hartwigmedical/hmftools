package com.hartwig.hmftools.lilac.evidence

import com.hartwig.hmftools.lilac.SequenceCount
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.ExpectedAlleles
import org.apache.logging.log4j.LogManager

class PhasedEvidenceFactory(private val minFragmentsPerAllele: Int, private val minFragmentsToRemoveSingles: Int, private val minEvidence: Int) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }


    fun evidence(context: HlaContext, fragments: List<AminoAcidFragment>): List<PhasedEvidence> {
        logger.info("Phasing HLA-${context.gene} records:")
        val result =  evidence(context.expectedAlleles, fragments)
        for (phasedEvidence in result) {
            logger.info("    $phasedEvidence")
        }

        return result
    }

    fun evidence(expectedAlleles: ExpectedAlleles, aminoAcidAminoAcidFragments: List<AminoAcidFragment>): List<PhasedEvidence> {
        val aminoAcidCounts = SequenceCount.aminoAcids(minEvidence, aminoAcidAminoAcidFragments)
        val heterozygousIndices = aminoAcidCounts.heterozygousLoci()
//        println(heterozygousIndices)

        val heterozygousEvidence = ExtendEvidence(minFragmentsPerAllele, minFragmentsToRemoveSingles, heterozygousIndices, aminoAcidAminoAcidFragments, expectedAlleles)

        val finalisedEvidence = mutableSetOf<PhasedEvidence>()
        val unprocessedEvidence = mutableListOf<PhasedEvidence>()
        unprocessedEvidence.addAll(heterozygousEvidence.pairedEvidence())


        while (unprocessedEvidence.isNotEmpty()) {
            val top = unprocessedEvidence.removeAt(0)
//            println("Processing: $top")

            val (parent, children) = heterozygousEvidence.merge(top, finalisedEvidence + unprocessedEvidence)

            if (children.isNotEmpty()) {
//                println("Produced:   $parent")
                finalisedEvidence.removeAll(children)
                unprocessedEvidence.removeAll(children)
                unprocessedEvidence.add(parent)
            } else {
                finalisedEvidence.add(parent)
            }

            unprocessedEvidence.sort()
        }

        return finalisedEvidence.sortedBy { it.aminoAcidIndices[0] }

    }

}


