package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair

class OncoSomaticEventReader {
    companion object {
        private val FUSION_SEPARATORS = setOf("-", " - ", "?")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ROS1", "CD74"),
                                            FusionPair("EP300", "MLL"),
                                            FusionPair("EP300", "MOZ"),
                                            FusionPair("RET", "CCDC6"))
        private val fusionReader = FusionReader(separators = FUSION_SEPARATORS, flipSet = FUSIONS_TO_FLIP)
    }

    fun read(gene: String, transcript: String, alteration: String): List<SomaticEvent> {
        val fusionsAndCnvs = listOfNotNull(readFusions(gene, alteration), readCnv(gene, alteration))
        return if (fusionsAndCnvs.isEmpty()) {
            listOfNotNull(ProteinAnnotation(transcript, alteration))
        } else {
            fusionsAndCnvs
        }
    }

    private fun readFusions(gene: String, alteration: String): FusionEvent? {
        return when {
            alteration.contains(Regex("Fusion")) -> fusionReader.read(gene, alteration)
            else                                 -> null
        }
    }

    private fun readCnv(gene: String, alteration: String): CnvEvent? {
        return when (alteration) {
            "Amplification" -> CnvEvent(gene, "Amplification")
            "Deletion"      -> CnvEvent(gene, "Deletion")
            else            -> null
        }
    }
}
