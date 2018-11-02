package com.hartwig.hmftools.knowledgebaseimporter.gene

import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.common.region.HmfTranscriptRegion
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericMutation
import org.apache.logging.log4j.LogManager

class GeneModel() {
    companion object {
        private val hmfGenePanel by lazy { HmfGenePanelSupplier.allGenesMap37() }
        private val logger = LogManager.getLogger("GeneModel")
    }

    fun hmfTranscriptForGene(gene: String): String? {
        // TODO (KODU): Deal with MLL2 (Can maybe map to KMT2B, KMT2A or something else) -> used by CGI
        val transcriptRegion = hmfGenePanel[gene]
        if (transcriptRegion == null) {
            logger.warn("Gene $gene not found in HMF gene panel!")
            return null
        }
        return transcriptRegion.transcriptID()
    }

    fun hmfTranscriptRegionForGenericMutation(mutation: GenericMutation): HmfTranscriptRegion? {
        // TODO (KODU): Deal with MLL2 (Can maybe map to KMT2B, KMT2A or something else) -> used by CGI
        val transcriptRegion = hmfGenePanel[mutation.gene]
        if (transcriptRegion == null) {
            logger.warn("Gene ${mutation.gene} not found in HMF gene panel!")
        } else if (mutation.transcript != null) {
            if (transcriptRegion.transcriptID() != mutation.transcript) {
                logger.warn("Non-canonical transcript ${mutation.transcript} used in knowledgebase for gene ${mutation.gene}")
                return null
            }
        }

        return transcriptRegion
    }
}
