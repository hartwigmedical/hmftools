package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.*
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarProteinAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import htsjdk.samtools.reference.IndexedFastaSequenceFile


class OncoKb(annotatedVariantsLocation: String, actionableVariantsLocation: String, transvarLocation: String,
             diseaseOntology: DiseaseOntology, private val reference: IndexedFastaSequenceFile) : Knowledgebase {
    companion object {
        private const val SOURCE: String = "oncoKb"
        private val FUSION_SEPARATORS = listOf("-", " - ", "?")
        private val FUSIONS_TO_FLIP = setOf(FusionPair("ROS1", "CD74"),
                                            FusionPair("EP300", "MLL"),
                                            FusionPair("EP300", "MOZ"),
                                            FusionPair("RET", "CCDC6"))
    }

    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)
    private val annotatedRecords by lazy { readTSVRecords(annotatedVariantsLocation) { OncoAnnotatedVariantRecord(it) }.map { preProcess(it) } }
    private val actionableRecords by lazy { readTSVRecords(actionableVariantsLocation) { OncoActionableVariantRecord(it) } }
    val cancerTypes by lazy { actionableRecords.map { it.cancerType }.map { Pair(it, diseaseOntology.findDoids(it)) }.toMap() }

    override val knownVariants: List<KnownVariantOutput> by lazy { knownVariants() }
    override val knownFusionPairs: List<FusionPair> by lazy { fusions().filterIsInstance<FusionPair>() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { fusions().filterIsInstance<PromiscuousGene>() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableVariants() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableCNVs() }
    override val actionableFusions: List<ActionableFusionOutput> by lazy { actionableFusions() }

    private fun knownVariants(): List<KnownVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(annotatedRecords.map { ProteinAnnotation(it.transcript, it.alteration) })
        return annotatedRecords.zip(transvarOutput)
                .flatMap { (oncoRecord, transvarOutput) ->
                    extractVariants(transvarOutput, reference)
                            .map {
                                KnownVariantOutput(oncoRecord.gene,
                                                   oncoRecord.transcript,
                                                   oncoRecord.oncogenicity,
                                                   SomaticVariantOutput(it))
                            }
                }
    }

    private fun fusions(): List<Fusion> {
        return annotatedRecords.filter { it.alteration.contains(Regex("Fusion")) }
                .map { extractFusion(it.gene, it.alteration, FUSION_SEPARATORS) }
                .map { flipFusion(it, FUSIONS_TO_FLIP) }
                .distinct()
    }

    private fun actionableFusions(): List<ActionableFusionOutput> {
        val fusionRecords = actionableRecords.filter { it.alteration.contains(Regex("Fusion")) }
        val fusions = fusionRecords.map { extractFusion(it.gene, it.alteration, FUSION_SEPARATORS) }.map { flipFusion(it, FUSIONS_TO_FLIP) }
        return fusionRecords.zip(fusions)
                .flatMap { (record, fusion) ->
                    record.drugs.map { drug ->
                        ActionableFusionOutput(fusion, actionability(drug, record))
                    }
                }
    }

    private fun actionableVariants(): List<ActionableVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(actionableRecords.map { ProteinAnnotation(it.transcript, it.alteration) })
        return actionableRecords.zip(transvarOutput)
                .flatMap { (record, transvarOutput) -> extractVariants(transvarOutput, reference).map { Pair(record, it) } }
                .flatMap { (record, somaticVariant) ->
                    record.drugs.map { drug ->
                        ActionableVariantOutput(record.gene, SomaticVariantOutput(somaticVariant), actionability(drug, record))
                    }
                }
    }

    private fun actionableCNVs(): List<ActionableCNVOutput> {
        return actionableRecords.filter { it.alteration == "Amplification" || it.alteration == "Deletion" }.flatMap { record ->
            record.drugs.map { drug ->
                ActionableCNVOutput(record.gene, record.alteration, actionability(drug, record))
            }
        }
    }

    private fun preProcess(record: OncoAnnotatedVariantRecord): OncoAnnotatedVariantRecord {
        return when {
            record.alteration.contains(Regex("IGH-NKX2")) && record.gene == "NKX2-1" ->
                record.copy(alteration = record.alteration.replace("IGH-NKX2", "IGH-NKX2-1"))
            else                                                                     -> record
        }
    }

    private fun actionability(drug: String, record: OncoActionableVariantRecord): Actionability {
        return Actionability(SOURCE, record.cancerType, drug, record.level, record.significance, "")
    }
}
