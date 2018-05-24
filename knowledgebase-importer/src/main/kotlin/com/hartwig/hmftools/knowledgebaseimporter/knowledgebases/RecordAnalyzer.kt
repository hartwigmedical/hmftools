package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.Annotation
import htsjdk.samtools.reference.IndexedFastaSequenceFile

class RecordAnalyzer(transvarLocation: String, private val reference: IndexedFastaSequenceFile) {
    private val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)

    fun knownVariants(knowledgebases: List<KnowledgebaseSource<*, *>>): List<KnownVariantOutput> {
        val records = knowledgebases.flatMap { it.knownKbRecords }
        val knownSomaticVariants = extractSomaticVariants(records)
        return knownSomaticVariants.map { (record, variant) ->
            KnownVariantOutput(record.transcript, record.additionalInfo, variant)
        }
    }

    fun actionableItems(knowledgebases: List<KnowledgebaseSource<*, *>>): List<ActionableKbItem> {
        val records = knowledgebases.flatMap { it.actionableKbRecords }
        val fusions = records.collectEvents<FusionEvent, ActionableRecord>()
        val cnvs = records.collectEvents<CnvEvent, ActionableRecord>()
        val somaticVariants = extractSomaticVariants(records)
        val events = fusions + cnvs + somaticVariants
        return events.flatMap { (record, event) ->
            record.actionability.map { ActionableKbItem(event, it) }
        }
    }

    private fun <R : KnowledgebaseRecord> extractSomaticVariants(records: List<R>): List<Pair<R, SomaticVariantEvent>> {
        return extractGdnaVariants(records) + extractTransvarVariants(records) + extractVariants(records)
    }

    private fun <R : KnowledgebaseRecord> extractGdnaVariants(records: List<R>): List<Pair<R, SomaticVariantEvent>> {
        val recordEventPairs = records.collectEvents<GDnaVariant, R>()
        return recordEventPairs.flatMap { (record, gdnaVariant) ->
            listOfNotNull(extractVariant(gdnaVariant.gDnaImpact, reference))
                    .map { Pair(record, SomaticVariantEvent(record.gene, it)) }
        }
    }

    private fun <R : KnowledgebaseRecord> extractTransvarVariants(records: List<R>): List<Pair<R, SomaticVariantEvent>> {
        return analyzeAnnotations(records, proteinAnalyzer) + analyzeAnnotations(records, cdnaAnalyzer)
    }

    private fun <R : KnowledgebaseRecord> extractVariants(records: List<R>): List<Pair<R, SomaticVariantEvent>> {
        return records.collectEvents<KnowledgebaseVariant, R>().map { Pair(it.first, extractVariant(it.second)) }
    }

    private fun extractVariant(variant: KnowledgebaseVariant): SomaticVariantEvent {
        val position = variant.position
        return when {
            variant.ref.isNullOrBlank() -> {
                val base = reference.getSubsequenceAt(variant.chromosome, position, position).baseString
                SomaticVariantEvent(variant.gene, variant.chromosome, position.toString(), base, base + variant.alt)
            }
            variant.alt.isNullOrBlank() -> {
                val base = reference.getSubsequenceAt(variant.chromosome, position - 1, position - 1).baseString
                SomaticVariantEvent(variant.gene, variant.chromosome, (position - 1).toString(), base + variant.ref, base)
            }
            else                        -> SomaticVariantEvent(variant.gene, variant.chromosome, position.toString(), variant.ref!!,
                                                               variant.alt!!)
        }
    }

    private inline fun <reified T, R : KnowledgebaseRecord> List<R>.collectEvents(): List<Pair<R, T>> {
        return this.flatMap { record -> record.events.filterIsInstance<T>().map { Pair(record, it) } }
    }

    private inline fun <reified T, R> analyzeAnnotations(records: List<R>, analyzer: TransvarAnalyzer<T>)
            : List<Pair<R, SomaticVariantEvent>> where T : Annotation, R : KnowledgebaseRecord {
        val recordEventPairs = records.collectEvents<T, R>()
        val transvarOutput = analyzer.analyze(recordEventPairs.map { it.second })
        return recordEventPairs.map { it.first }.zip(transvarOutput).flatMap { (record, output) ->
            val somaticVariants = extractVariants(output, reference)
            somaticVariants.map { Pair(record, SomaticVariantEvent(record.gene, it)) }
        }
    }
}
