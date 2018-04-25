package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
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
        private const val FUSION_SEPARATOR = "-"
    }

    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)
    private val annotatedRecords = readCSVRecords(annotatedVariantsLocation) { OncoAnnotatedVariantRecord(it) }
    private val actionableRecords = readCSVRecords(actionableVariantsLocation) { OncoActionableVariantRecord(it) }

    override val knownVariants: List<KnownVariantOutput> = knownVariants()
    override val knownFusionPairs: List<Pair<String, String>>
        get() = TODO("not implemented") //To change initializer of created properties use File | Settings | File Templates.
    override val promiscuousGenes: List<String>
        get() = TODO("not implemented") //To change initializer of created properties use File | Settings | File Templates.
    override val actionableVariants: List<ActionableVariantOutput> = actionableVariants()
    override val actionableCNVs: List<ActionableCNVOutput> = actionableCNVs()
    override val actionableFusions: List<ActionableFusionOutput>
        get() = TODO("not implemented") //To change initializer of created properties use File | Settings | File Templates.

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

    private fun knownVariants(): List<KnownVariantOutput> {
        val transvarOutput = proteinAnalyzer.analyze(annotatedRecords.map { ProteinAnnotation(it.transcript, it.alteration) })
        return annotatedRecords.zip(transvarOutput)
                .flatMap { (oncoRecord, transvarOutput) ->
                    extractVariants(transvarOutput, reference)
                            .map { KnownVariantOutput(oncoRecord.gene, oncoRecord.transcript, oncoRecord.oncogenicity, it) }
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
}
