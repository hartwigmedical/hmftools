package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ExonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericMutation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.GeneMutations

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
            isFusion(alteration)          -> listOfNotNull(readFusions(gene, alteration, effect))
            isCnv(alteration)             -> listOfNotNull(readCnv(gene, alteration))
            isGenericMutation(alteration) -> listOfNotNull(readGenericMutation(gene, transcript, alteration))
            else                          -> listOf(ProteinAnnotation(transcript, alteration))
        }
    }

    private fun isFusion(alteration: String): Boolean = alteration.contains("Fusion".toRegex())
    private fun isCnv(alteration: String): Boolean = alteration == "Amplification" || alteration == "Deletion"
    private fun isGenericMutation(alteration: String): Boolean = alteration == "Oncogenic Mutations" || isExonMutation(alteration)
    private fun isExonMutation(alteration: String) = alteration.contains("${EXON_PATTERN}Mutations".toRegex(RegexOption.IGNORE_CASE))

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
