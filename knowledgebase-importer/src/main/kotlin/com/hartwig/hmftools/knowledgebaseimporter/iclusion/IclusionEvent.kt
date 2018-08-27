package com.hartwig.hmftools.knowledgebaseimporter.iclusion

import com.hartwig.hmftools.apiclients.iclusion.data.IclusionMutationDetails
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers.IclusionExonMutationReader.EXON_MUTATION_REGEX
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers.IclusionFusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.EventType
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.TransvarMatcher
import org.apache.logging.log4j.LogManager

data class IclusionEvent(override val gene: String, override val variant: String, val transcript: String,
                         override val types: List<EventType> = emptyList()) : KnowledgebaseEvent {
    companion object {
        private val logger = LogManager.getLogger("IclusionEvent")

        operator fun invoke(mutation: IclusionMutationDetails, transcript: String): IclusionEvent {
            return IclusionEvent(mutation.geneName, mutation.variantName, transcript, eventTypes(mutation.geneName, mutation.variantName))
        }

        private fun eventTypes(gene: String, variant: String): List<EventType> {
            val types = mutableSetOf<EventType>()
            types.addAll(cnvTypes(variant))
            types.addAll(mutationTypes(variant))
            types.addAll(expressionTypes(variant))
            types.addAll(fusionTypes(gene, variant))
            if (variant == "Wild-type") types.add(EventType.WILD_TYPE)
            if (gene == "MSI") types.add(EventType.MSI)
            if (types.isEmpty()) {
                logger.warn("could not map: $gene $variant to event type")
            }
            return types.toList()
        }

        private fun cnvTypes(variant: String): List<EventType> {
            return if (variant == "COPY-GAIN" || variant == "DELETION" || variant == "AMPLIFICATION" || variant == "LOSS" ||
                    variant == "COPY NUMBER VARIATION") {
                listOf(EventType.CNV)
            } else {
                emptyList()
            }
        }

        private fun mutationTypes(variant: String): List<EventType> {
            return if (variant == "MUTATION" || variant == "ALTERATION" || variant == "ACTIVATING MUTATION" || variant == "INACTIVATING MUTATION" ||
                    variant.contains("Splicing alteration") || TransvarMatcher.contains(variant) || variant.matches(EXON_MUTATION_REGEX) ||
                    variant.toLowerCase().contains("utr mutation") || variant.toLowerCase().contains("utr alteration") || variant == "TRUNCATING MUTATION") {
                listOf(EventType.MUT)
            } else {
                emptyList()
            }
        }

        private fun expressionTypes(variant: String): List<EventType> {
            return if (variant == "OVEREXPRESSION" || variant == "UNDEREXPRESSION" || variant == "EXPRESSION") {
                listOf(EventType.EXPR)
            } else {
                emptyList()
            }
        }

        private fun fusionTypes(gene: String, variant: String): List<EventType> {
            return if (variant.toLowerCase() == "fusions" || variant.toLowerCase().contains("fusion") ||
                    (variant.contains(IclusionFusionReader.FUSION_SEPARATOR) && gene.length > 2 &&
                            variant.toLowerCase().contains(gene.toLowerCase().substring(0, 3)))) {
                return listOf(EventType.FUS)
            } else {
                emptyList()
            }
        }
    }
}
