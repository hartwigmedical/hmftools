package com.hartwig.hmftools.knowledgebaseimporter.cosmic

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.Doid
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import org.apache.commons.csv.CSVRecord

class Cosmic(fusionsLocation: String) : Knowledgebase {
    override val source = "cosmic"
    override val knownVariants: List<KnownVariantOutput> = listOf()
    override val knownFusionPairs by lazy { readCSVRecords(fusionsLocation) { readAndCorrectFusion(it) }.filterNotNull().distinct() }
    override val promiscuousGenes: List<PromiscuousGene> = listOf()
    override val actionableVariants: List<ActionableVariantOutput> = listOf()
    override val actionableCNVs: List<ActionableCNVOutput> = listOf()
    override val actionableFusionPairs: List<ActionableFusionPairOutput> = listOf()
    override val actionablePromiscuousGenes: List<ActionablePromiscuousGeneOutput> = listOf()
    override val actionableRanges: List<ActionableGenomicRangeOutput> = listOf()
    override val cancerTypes: Map<String, Set<Doid>> = mapOf()

    private fun readAndCorrectFusion(csvRecord: CSVRecord): FusionPair? {
        val fiveGene = curateGene(csvRecord["5' Partner"].split("_").first())
        val threeGene = curateGene(csvRecord["3' Partner"].split("_").first())

        if (!isBlackListed(fiveGene, threeGene)) {
            return FusionPair(fiveGene, threeGene)
        }

        return null
    }

    private fun curateGene(gene: String): String {
        return when (gene) {
            "DUX4L1" -> "DUX4"
            "SIP1" -> "GEMIN2"
            "KIAA0284" -> "CEP170B"
            "ACCN1" -> "ASIC2"
            "FAM22A" -> "NUTM2A"
            "ROD1" -> "PTBP3"
            else -> gene
        }
    }

    private fun isBlackListed(fiveGene: String, threeGene: String) : Boolean {
        // The below 6 fusions all entered COSMIC based on a single publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3398135/
        if (fiveGene == "INTS4" && threeGene == "GAB2") {
            return true
        } else if (fiveGene == "PLA2R1" && threeGene == "RBMS1") {
            return true
        } else if (fiveGene == "FBXL18" && threeGene == "RNF216") {
            return true
        } else if (fiveGene == "PLXND1" && threeGene == "TMCC1") {
            return true
        } else if (fiveGene == "SEPT8" && threeGene == "AFF4") {
            return true
        } else if (fiveGene == "IL6R" && threeGene == "ATP8B2") {
            return true
        }

        // See DEV-1061. Found once in a single publication: https://www.ncbi.nlm.nih.gov/pubmed/20033038
        if (fiveGene == "NF1" && threeGene == "ASIC2") {
            return true
        }

        return false
    }
}
