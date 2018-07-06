package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair

class OncoSomaticEventReader {
    companion object {
        private const val EXON_PATTERN = "Exon\\s?([0-9]+)\\s+"
        private val FUSION_SEPARATORS = setOf("-", " - ", "?")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ROS1", "CD74"),
                                            FusionPair("EP300", "MLL"),
                                            FusionPair("EP300", "MOZ"),
                                            FusionPair("RET", "CCDC6"))
        private val fusionReader = FusionReader(separators = FUSION_SEPARATORS, flipSet = FUSIONS_TO_FLIP)
    }

    fun read(gene: String, transcript: String, alteration: String, effect: String = ""): List<SomaticEvent> {
        return when {
            isFusion(alteration)                                      -> listOfNotNull(readFusions(gene, alteration, effect))
            isCnv(alteration)                                         -> listOfNotNull(readCnv(gene, alteration))
            isTruncatingMutation(alteration)                          -> emptyList()
            isMSI(alteration)                                         -> emptyList()
            isWildtype(alteration)                                    -> emptyList()
            isSplice(alteration)                                      -> emptyList()
            isExonInsertion(alteration) && isExonDeletion(alteration) -> emptyList()
            isExonInsertion(alteration)                               -> emptyList()
            isExonDeletion(alteration)                                -> emptyList()
            isGainOfFunction(alteration)                              -> emptyList()
            isInternalTandemDuplication(alteration)                   -> emptyList()
            isBrafV600EV600K(alteration)                              -> listOf(ProteinAnnotation(transcript, "V600E"),
                                                                                ProteinAnnotation(transcript, "V600K"))
            isGenericMutation(alteration)                             -> listOfNotNull(readGenericMutation(gene, transcript, alteration))
            else                                                      -> listOf(ProteinAnnotation(transcript, alteration))
        }
    }

    private fun isFusion(alteration: String): Boolean = alteration.contains("Fusion".toRegex())
    private fun isCnv(alteration: String): Boolean = alteration == "Amplification" || alteration == "Deletion"
    private fun isGenericMutation(alteration: String): Boolean = alteration == "Oncogenic Mutations" || isExonMutation(alteration)
    private fun isExonMutation(alteration: String) = alteration.contains("${EXON_PATTERN}Mutations".toRegex(RegexOption.IGNORE_CASE))
    private fun isTruncatingMutation(alteration: String) = alteration == "Truncating Mutations"
    private fun isMSI(alteration: String) = alteration == "Microsatellite Instability-High"
    private fun isWildtype(alteration: String) = alteration == "Wildtype"
    private fun isSplice(alteration: String) = alteration.contains("splice")
    private fun isGainOfFunction(alteration: String) = alteration == "Gain-of-function Mutations"
    private fun isBrafV600EV600K(alteration: String) = alteration == "p61BRAF-V600E"
    private fun isInternalTandemDuplication(alteration: String) = alteration.contains("internal tandem duplications")
    private fun isExonInsertion(alteration: String) =
            alteration.contains(EXON_PATTERN.toRegex(RegexOption.IGNORE_CASE)) && alteration.contains("insertion")

    private fun isExonDeletion(alteration: String) =
            alteration.contains(EXON_PATTERN.toRegex(RegexOption.IGNORE_CASE)) && alteration.contains("deletion")


    private fun readFusions(gene: String, alteration: String, effect: String): FusionEvent? {
        return when {
            effect.contains("Loss-of-function") -> null
            else                                -> fusionReader.read(gene, alteration)
        }
    }

    private fun readCnv(gene: String, alteration: String): CnvEvent? {
        return when (alteration) {
            "Amplification" -> CnvEvent(gene, "Amplification")
            "Deletion"      -> CnvEvent(gene, "Deletion")
            else            -> null
        }
    }

    private fun readGenericMutation(gene: String, transcript: String, alteration: String): GenericMutation? {
        return when {
            alteration == "Oncogenic Mutations" -> GeneMutations(gene, transcript)
            isExonMutation(alteration)          -> ExonMutations(gene, transcript, extractExon(alteration))
            else                                -> null
        }
    }

    private fun extractExon(alteration: String): Int {
        val matchResult = EXON_PATTERN.toRegex(RegexOption.IGNORE_CASE).find(alteration)
        return matchResult!!.groupValues[1].toInt()
    }
}
