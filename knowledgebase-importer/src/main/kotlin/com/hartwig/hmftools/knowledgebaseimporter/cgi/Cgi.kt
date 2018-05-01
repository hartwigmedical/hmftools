package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.*
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import htsjdk.samtools.reference.IndexedFastaSequenceFile

class Cgi(variantsLocation: String, biomarkersLocation: String, transvarLocation: String, private val reference: IndexedFastaSequenceFile) :
        Knowledgebase {
    companion object {
        private const val SOURCE: String = "cgi"
        private val FUSION_SEPARATORS = listOf("__")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ABL1", "BCR"),
                                            FusionPair("PDGFRA", "FIP1L1"),
                                            FusionPair("PDGFB", "COL1A1"))
        private val FUSIONS_TO_FILTER = setOf(FusionPair("RET", "TPCN1"))
    }

    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)
    private val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    private val somaticVariantRecords = readCSVRecords(variantsLocation) { CgiKnownVariantRecord(it) }.filter { it.context == "somatic" }
    private val biomarkersRecords = readCSVRecords(biomarkersLocation) { CgiBiomarkersRecord(it) }

    override val knownVariants: List<KnownVariantOutput> by lazy { knownVariants() }
    override val knownFusionPairs: List<FusionPair> by lazy { fusionRecords().filterIsInstance<FusionPair>() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { fusionRecords().filterIsInstance<PromiscuousGene>() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableVariants() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableCNVs() }
    override val actionableFusions: List<ActionableFusionOutput>
        get() = TODO("not implemented") //To change initializer of created properties use File | Settings | File Templates.

    private fun knownVariants(): List<KnownVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(somaticVariantRecords.map { ProteinAnnotation(it.transcript, it.impact) })
        return somaticVariantRecords.zip(transvarOutput)
                .flatMap { (cgiRecord, transvarOutput) ->
                    val cgiVariants = extractCgiVariants(cgiRecord.gdna, reference)
                    val inferredKnownVariants = extractVariants(transvarOutput, reference)
                            .filterNot { variant -> cgiVariants.any { it == variant } }
                            .map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "", it) }
                    val cgiKnownVariants = cgiVariants.map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "CGI", it) }
                    cgiKnownVariants + inferredKnownVariants
                }
    }

    private fun actionableVariants(): List<ActionableVariantOutput> {
        val transvarOutput = annotateCgiBiomarkers()
        return biomarkersRecords.zip(transvarOutput)
                .flatMap { (cgiRecord, transvarOutput) ->
                    val cgiVariants = extractCgiVariants(cgiRecord.gdna, reference)
                    val inferredProteinVariants = extractVariants(transvarOutput.first, reference)
                            .filterNot { variant -> cgiVariants.any { it == variant } }
                    val inferredCdnaVariants = extractVariants(transvarOutput.second, reference)
                            .filterNot { variant -> cgiVariants.any { it == variant } || inferredProteinVariants.any { it == variant } }
                    actionableVariantOutput(cgiRecord, cgiVariants + inferredProteinVariants + inferredCdnaVariants)
                }
    }

    private fun actionableCNVs(): List<ActionableCNVOutput> {
        val cnaRecords = biomarkersRecords.filter { it.alterationType == "CNA" }
        return cnaRecords.flatMap { record ->
            record.cancerTypes.map { cancerType ->
                ActionableCNVOutput(record.gene, extractAmpOrDel(record.alteration), actionability(cancerType, record))
            }
        }
    }

    private fun fusionRecords(): List<Fusion> {
        return biomarkersRecords.filter { it.alterationType == "FUS" }
                .map { extractFusion(it.gene, it.alteration.trim(), FUSION_SEPARATORS) }
                .map { flipFusion(it, FUSIONS_TO_FLIP) }
                .filterNot { FUSIONS_TO_FILTER.contains(it) }
                .distinct()
    }

    private fun extractCgiVariants(gdna: String, reference: IndexedFastaSequenceFile): List<SomaticVariant> {
        val cgiVariantGdnas = gdna.split("__").filterNot { it.isBlank() }
        val chromosome = extractChromosome(cgiVariantGdnas.first())
        return extract(chromosome, cgiVariantGdnas, reference)
    }

    private fun annotateCgiBiomarkers(): List<Pair<TransvarOutput, TransvarOutput>> {
        val proteinOutput = proteinAnalyzer.analyze(biomarkersRecords.map { ProteinAnnotation(it.transcript, it.protein) })
        val cdnaOutput = cdnaAnalyzer.analyze(biomarkersRecords.map { CDnaAnnotation(it.transcript, it.cdna) })
        return proteinOutput.zip(cdnaOutput)
    }

    private fun actionableVariantOutput(cgiRecord: CgiBiomarkersRecord,
                                        somaticVariants: List<SomaticVariant>): List<ActionableVariantOutput> {
        return cgiRecord.cancerTypes.flatMap { cancerType ->
            somaticVariants.map { somaticVariant ->
                ActionableVariantOutput(cgiRecord.gene, somaticVariant, actionability(cancerType, cgiRecord))
            }
        }
    }

    private fun extractAmpOrDel(alteration: String): String = if (alteration.contains("amp")) "Amplification" else "Deletion"

    private fun actionability(cancerType: String, cgiBiomarker: CgiBiomarkersRecord): Actionability {
        return Actionability(SOURCE, cancerType, cgiBiomarker.drug, cgiBiomarker.level, cgiBiomarker.association, "")
    }
}
