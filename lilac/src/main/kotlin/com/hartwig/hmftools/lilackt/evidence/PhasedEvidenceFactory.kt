package com.hartwig.hmftools.lilackt.evidence

import com.hartwig.hmftools.lilackt.LilacConfig
import com.hartwig.hmftools.lilackt.SequenceCount
import com.hartwig.hmftools.lilackt.amino.AminoAcidFragment
import com.hartwig.hmftools.lilackt.hla.HlaContext
import com.hartwig.hmftools.lilackt.nuc.ExpectedAlleles
import org.apache.logging.log4j.LogManager

class PhasedEvidenceFactory(private val config: LilacConfig) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }


    fun evidence(context: HlaContext, fragments: List<AminoAcidFragment>): List<PhasedEvidence> {
        logger.info("Phasing HLA-${context.gene} records:")
        val result =  evidence(context.expectedAlleles, fragments)

        if (config.debugPhasing) {
            logger.info("    Consolidating evidence")
        }

        for (phasedEvidence in result) {
            logger.info("    $phasedEvidence")
        }

        return result
    }

    fun evidence(expectedAlleles: ExpectedAlleles, aminoAcidAminoAcidFragments: List<AminoAcidFragment>): List<PhasedEvidence> {
        val aminoAcidCounts = SequenceCount.aminoAcids(config.minEvidence, aminoAcidAminoAcidFragments)
        val heterozygousIndices = aminoAcidCounts.heterozygousLoci()
        if (config.debugPhasing) {
            logger.info("    Heterozygous Indices: $heterozygousIndices")
        }

        val heterozygousEvidence = ExtendEvidence(config, heterozygousIndices, aminoAcidAminoAcidFragments, expectedAlleles)

        val finalisedEvidence = mutableSetOf<PhasedEvidence>()
        val unprocessedEvidence = mutableListOf<PhasedEvidence>()
        unprocessedEvidence.addAll(heterozygousEvidence.pairedEvidence())

        if (config.debugPhasing) {
            logger.info("    Extending paired evidence")
        }


        while (unprocessedEvidence.isNotEmpty()) {
            val top = unprocessedEvidence.removeAt(0)
            if (config.debugPhasing) {
                logger.info("    Processing top: $top")
            }

            val (parent, children) = heterozygousEvidence.merge(top, finalisedEvidence + unprocessedEvidence)

            if (children.isNotEmpty()) {
                if (config.debugPhasing) {
                    logger.info("    Produced child: $parent")
                }

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


