package com.hartwig.hmftools.knowledgebaseimporter.cosmic

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.Doid
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import org.apache.commons.csv.CSVRecord

class Cosmic(fusionsLocation: String) : Knowledgebase {
    override val source = "cosmic"
    override val knownVariants: List<KnownVariantOutput> = listOf()
    override val knownFusionPairs by lazy { readCSVRecords(fusionsLocation) { readAndCorrectFusion(it) }.distinct() }
    override val promiscuousGenes: List<PromiscuousGene> = listOf()
    override val actionableVariants: List<ActionableVariantOutput> = listOf()
    override val actionableCNVs: List<ActionableCNVOutput> = listOf()
    override val actionableFusionPairs: List<ActionableFusionPairOutput> = listOf()
    override val actionablePromiscuousGenes: List<ActionablePromiscuousGeneOutput> = listOf()
    override val actionableRanges: List<ActionableGenomicRangeOutput> = listOf()
    override val cancerTypes: Map<String, Set<Doid>> = mapOf()

    private fun readAndCorrectFusion(csvRecord: CSVRecord): FusionPair {
        val fiveGene = csvRecord["5' Partner"].split("_").first()
        val threeGene = csvRecord["3' Partner"].split("_").first()
        if (threeGene.equals("DUX4L1")) {
            threeGene.replace("DUX4L1", "DUX4")
        }
        if (threeGene.equals("SIP1")) {
            threeGene.replace("SIP1", "GEMIN2")
        }
        if (threeGene.equals("KIAA0284")) {
            threeGene.replace("KIAA0284", "CEP170B")
        }
        if (threeGene.equals("ACCN1")) {
            threeGene.replace("ACCN1", "ASIC2")
        }
        if (threeGene.equals("FAM22A")) {
            threeGene.replace("FAM22A", "NUTM2A")
        }
        return FusionPair(fiveGene, threeGene)
    }
}
