package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.extractFusion
import com.hartwig.hmftools.knowledgebaseimporter.flipFusion
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation

class OncoSomaticEventReader {
    companion object {
        private val FUSION_SEPARATORS = listOf("-", " - ", "?")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ROS1", "CD74"),
                                            FusionPair("EP300", "MLL"),
                                            FusionPair("EP300", "MOZ"),
                                            FusionPair("RET", "CCDC6"))
        private val FUSIONS_TO_FILTER = setOf<FusionPair>()
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
        val fusion = when {
            alteration.contains(Regex("Fusion")) -> extractFusion(gene, alteration, FUSION_SEPARATORS)
            else                                 -> null
        }
        return if (FUSIONS_TO_FILTER.contains(fusion) || fusion == null) null else flipFusion(fusion, FUSIONS_TO_FLIP)
    }

    private fun readCnv(gene: String, alteration: String): CnvEvent? {
        return when (alteration) {
            "Amplification" -> CnvEvent(gene, "Amplification")
            "Deletion"      -> CnvEvent(gene, "Deletion")
            else            -> null
        }
    }
}
