package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.extractFusion
import com.hartwig.hmftools.knowledgebaseimporter.flipFusion
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readTSVRecords
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import htsjdk.samtools.reference.IndexedFastaSequenceFile

class Cgi(variantsLocation: String, biomarkersLocation: String, transvarLocation: String, diseaseOntology: DiseaseOntology,
          private val reference: IndexedFastaSequenceFile) :
        Knowledgebase {
    companion object {
        private val FUSION_SEPARATORS = listOf("__")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ABL1", "BCR"),
                                            FusionPair("PDGFRA", "FIP1L1"),
                                            FusionPair("PDGFB", "COL1A1"))
        private val FUSIONS_TO_FILTER = setOf(FusionPair("RET", "TPCN1"))
    }

    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)
    private val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    private val somaticVariantRecords by lazy { readTSVRecords(variantsLocation) { CgiKnownVariantRecord(it) }.filter { it.context == "somatic" } }
    private val biomarkersRecords by lazy { readTSVRecords(biomarkersLocation) { CgiBiomarkersRecord(it) } }
    val cancerTypes by lazy { biomarkersRecords.flatMap { it.cancerTypes }.map { Pair(it, diseaseOntology.findDoids(it)) }.toMap() }

    override val source = "cgi"
    override val knownVariants: List<KnownVariantOutput> by lazy { knownVariants() }
    override val knownFusionPairs: List<FusionPair> by lazy { actionableFusions.map { it.fusion }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { actionableFusions.map { it.fusion }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableVariants() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableCNVs() }
    override val actionableFusions: List<ActionableFusionOutput> by lazy { actionableFusions() }

    private fun knownVariants(): List<KnownVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(somaticVariantRecords.map { ProteinAnnotation(it.transcript, it.impact) })
        return somaticVariantRecords.zip(transvarOutput)
                .flatMap { (cgiRecord, transvarOutput) ->
                    val cgiVariants = extractCgiVariants(cgiRecord.gdna, reference)
                    val inferredKnownVariants = extractVariants(transvarOutput, reference)
                            .filterNot { variant -> cgiVariants.any { it == variant } }
                            .map { KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "", SomaticVariantOutput(it)) }
                    val cgiKnownVariants = cgiVariants.map {
                        KnownVariantOutput(cgiRecord.gene, cgiRecord.transcript, "CGI", SomaticVariantOutput(it))
                    }
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

    private fun actionableFusions(): List<ActionableFusionOutput> {
        val fusionRecords = biomarkersRecords.filter { it.alterationType == "FUS" }
        val fusions = fusionRecords.map { extractFusion(it.gene, it.alteration.trim(), FUSION_SEPARATORS) }
                .map { flipFusion(it, FUSIONS_TO_FLIP) }
        return fusionRecords.zip(fusions).filterNot { FUSIONS_TO_FILTER.contains(it.second) }
                .flatMap { (record, fusion) ->
                    record.cancerTypes.map { cancerType -> ActionableFusionOutput(fusion, actionability(cancerType, record)) }
                }
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
            somaticVariants.map {
                ActionableVariantOutput(cgiRecord.gene, SomaticVariantOutput(it), actionability(cancerType, cgiRecord))
            }
        }
    }

    private fun extractAmpOrDel(alteration: String): String = if (alteration.contains("amp")) "Amplification" else "Deletion"

    private fun actionability(cancerType: String, cgiBiomarker: CgiBiomarkersRecord): Actionability {
        return Actionability(source, cancerType, cgiBiomarker.drug, cgiBiomarker.level, cgiBiomarker.association, "")
    }
}
