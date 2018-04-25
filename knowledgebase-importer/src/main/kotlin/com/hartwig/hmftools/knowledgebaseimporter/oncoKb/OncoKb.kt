package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.extractFusion
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarProteinAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import htsjdk.samtools.reference.IndexedFastaSequenceFile


class OncoKb(annotatedVariantsLocation: String, actionableVariantsLocation: String,
             transvarLocation: String, private val reference: IndexedFastaSequenceFile) : Knowledgebase {
    companion object {
        private const val SOURCE: String = "oncoKb"
        private val FUSION_SEPARATORS = listOf("-", " - ", "?")
    }

    private val proteinAnalyzer by lazy { TransvarProteinAnalyzer(transvarLocation) }
    private val annotatedRecords by lazy { readCSVRecords(annotatedVariantsLocation) { OncoAnnotatedVariantRecord(it) }.map { preProcess(it) } }
    private val actionableRecords by lazy { readCSVRecords(actionableVariantsLocation) { OncoActionableVariantRecord(it) } }

    override val knownVariants: List<KnownVariantOutput> by lazy { knownVariants() }
    override val knownFusionPairs: List<Pair<String, String>> by lazy { knownFusions() }
    override val promiscuousGenes: List<String>
        get() = TODO("not implemented") //To change initializer of created properties use File | Settings | File Templates.
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableVariants() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableCNVs() }
    override val actionableFusions: List<ActionableFusionOutput>
        get() = TODO("not implemented") //To change initializer of created properties use File | Settings | File Templates.

    private fun knownVariants(): List<KnownVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(annotatedRecords.map { ProteinAnnotation(it.transcript, it.alteration) })
        return annotatedRecords.zip(transvarOutput)
                .flatMap { (oncoRecord, transvarOutput) ->
                    extractVariants(transvarOutput, reference)
                            .map { KnownVariantOutput(oncoRecord.gene, oncoRecord.transcript, oncoRecord.oncogenicity, it) }
                }
    }

    private fun knownFusions(): List<Pair<String, String>> {
        return annotatedRecords.filter { it.alteration.contains(Regex("Fusion$")) }.mapNotNull {
            extractFusion(it.gene, it.alteration.substringBefore("Fusion").trim(), FUSION_SEPARATORS)
        }.map { flipGenePairs(it) }.distinct()
    }

    private fun flipGenePairs(fusionPair: Pair<String, String>): Pair<String, String> {
        println(fusionPair)
        return when (fusionPair) {
            Pair("ROS1", "CD74") -> Pair("CD74", "ROS1")
            Pair("EP300", "MLL") -> Pair("MLL", "EP300")
            Pair("EP300", "MOZ") -> Pair("MOZ", "EP300")
            Pair("RET", "CCDC6") -> Pair("CCDC6", "RET")
            else                 -> fusionPair
        }
    }

    private fun actionableVariants(): List<ActionableVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(actionableRecords.map { ProteinAnnotation(it.transcript, it.alteration) })
        return actionableRecords.zip(transvarOutput)
                .flatMap { (record, transvarOutput) -> extractVariants(transvarOutput, reference).map { Pair(record, it) } }
                .flatMap { (record, somaticVariant) ->
                    record.drugs.map { drug ->
                        ActionableVariantOutput(record.gene,
                                                somaticVariant,
                                                Actionability(SOURCE, record.cancerType, drug, record.level, record.significance, ""))
                    }
                }
    }

    private fun actionableCNVs(): List<ActionableCNVOutput> {
        return actionableRecords.filter { it.alteration == "Amplification" || it.alteration == "Deletion" }.flatMap { record ->
            record.drugs.map { drug ->
                ActionableCNVOutput(record.gene,
                                    record.alteration,
                                    Actionability(SOURCE, record.cancerType, drug, record.level, record.significance, ""))
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

}
