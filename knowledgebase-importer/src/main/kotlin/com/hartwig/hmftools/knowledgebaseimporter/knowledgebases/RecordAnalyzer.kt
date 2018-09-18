package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import com.hartwig.hmftools.common.region.HmfTranscriptRegion
import com.hartwig.hmftools.knowledgebaseimporter.gene.GeneModel
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.GenomicRangeEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import htsjdk.samtools.reference.IndexedFastaSequenceFile

class RecordAnalyzer(transvarLocation: String, private val reference: IndexedFastaSequenceFile, private val geneModel: GeneModel) {
    companion object {
        private val blacklistedDrugs = setOf("chemotherapy", "aspirin", "steroid")
    }

    private val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)

    fun knownVariants(knowledgebases: List<KnowledgebaseSource<*, *>>): List<KnownVariantOutput> {
        val records = knowledgebases.flatMap { it.knownKbRecords }
        val somaticVariantRecords = extractSomaticVariants(records) + extractOtherSomaticVariants(records)
        val knownSomaticVariants = somaticVariantRecords.filter { it.second is SomaticVariantEvent }
                .map { Pair(it.first, it.second as SomaticVariantEvent) }
        return knownSomaticVariants.map { (record, variant) ->
            KnownVariantOutput(record.transcript, record.additionalInfo, variant)
        }
    }

    fun actionableItems(knowledgebases: List<KnowledgebaseSource<*, *>>): List<ActionableItem<*>> {
        val records = knowledgebases.flatMap { it.actionableKbRecords }
        val genomicRangeEvents = analyzeGenomicRanges(records)
        val somaticVariants = extractSomaticVariants(records)
        val actionableEvents = records.collectEvents<ActionableEvent, ActionableRecord>()
        val events = actionableEvents + somaticVariants + genomicRangeEvents
        return events.flatMap { (record, event) ->
            record.actionability.filter { !blacklistedDrugs.contains(it.drug.name.toLowerCase()) }.map { ActionableItem(event, it) }
        }
    }

    private fun <R : KnowledgebaseRecord> extractOtherSomaticVariants(records: List<R>): List<Pair<R, ActionableEvent>> {
        return extractGdnaVariants(records.collectOtherEvents()) + analyzeAnnotations(records.collectOtherEvents(), proteinAnalyzer) +
                analyzeAnnotations(records.collectEvents(), cdnaAnalyzer) + extractVariants(records.collectOtherEvents())
    }

    private fun <R : KnowledgebaseRecord> extractSomaticVariants(records: List<R>): List<Pair<R, ActionableEvent>> {
        return extractGdnaVariants(records.collectEvents()) + extractTransvarVariants(records) + extractVariants(records.collectEvents())
    }

    private fun <R : KnowledgebaseRecord> extractGdnaVariants(
            recordEventPairs: List<Pair<R, GDnaVariant>>): List<Pair<R, ActionableEvent>> {
        return recordEventPairs.flatMap { (record, gdnaVariant) ->
            listOfNotNull(extractVariant(record.gene, record.transcript, gdnaVariant.gDnaImpact, reference))
                    .map { Pair(record, it) }
        }
    }

    private fun <R : KnowledgebaseRecord> extractTransvarVariants(records: List<R>): List<Pair<R, ActionableEvent>> {
        return analyzeAnnotations(records.collectEvents(), proteinAnalyzer) + analyzeAnnotations(records.collectEvents(), cdnaAnalyzer)
    }

    private fun <R : KnowledgebaseRecord> extractVariants(
            recordEventPairs: List<Pair<R, KnowledgebaseVariant>>): List<Pair<R, SomaticVariantEvent>> {
        return recordEventPairs.map { Pair(it.first, extractVariant(it.second)) }
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
            else -> SomaticVariantEvent(variant.gene, variant.chromosome, position.toString(), variant.ref!!,
                    variant.alt!!)
        }
    }

    private inline fun <reified T, R : KnowledgebaseRecord> List<R>.collectOtherEvents(): List<Pair<R, T>> {
        return this.flatMap { record ->
            record.events.filterIsInstance<OtherEvents>().flatMap { it.events.filterIsInstance<T>() }.map { Pair(record, it) }
        }
    }

    private inline fun <reified T, R : KnowledgebaseRecord> List<R>.collectEvents(): List<Pair<R, T>> {
        return this.flatMap { record -> record.events.filterIsInstance<T>().map { Pair(record, it) } }
    }

    private inline fun <reified T, R> analyzeAnnotations(recordEventPairs: List<Pair<R, T>>, analyzer: TransvarAnalyzer<T>)
            : List<Pair<R, ActionableEvent>> where T : HgvsAnnotation, R : KnowledgebaseRecord {
        val transvarOutput = analyzer.analyze(recordEventPairs.map { it.second })
        return recordEventPairs.map { it.first }.zip(transvarOutput).flatMap { (record, output) ->
            extractVariants(record.gene, record.transcript, output, reference).map { Pair(record, it) }
        }
    }

    private fun <R : ActionableRecord> analyzeGenomicRanges(records: List<R>): List<Pair<R, GenomicRangeEvent>> {
        val genericMutationRecords = records.collectEvents<GenericMutation, R>()
        val geneModel = createGeneModel(genericMutationRecords.map { it.second })
        return genericMutationRecords.flatMap { (record, mutation) ->
            val gene = geneModel[mutation.gene]
            mutationCodingRange(mutation, gene).map {
                Pair(record, GenomicRangeEvent(mutation.gene, mutation.transcript.orEmpty(), gene!!.chromosome(), it.start.toString(),
                        it.endInclusive.toString(), gene.transcriptID()))
            }
        }
    }

    private fun mutationCodingRange(mutation: GenericMutation, gene: HmfTranscriptRegion?): List<ClosedRange<Long>> {
        if (gene == null) {
            return emptyList()
        }
        return when (mutation) {
            is GeneMutations -> listOf(gene.start()..gene.end())
            is ExonMutations -> listOf(gene.exome()[mutation.exonNumber].start()..gene.exome()[mutation.exonNumber].end())
            is CodonRangeMutations -> gene.codonRangeByIndex(mutation.startCodon, mutation.endCodon)?.map { it -> it.start()..it.end() }
                    ?: emptyList()
            is CodonMutations -> gene.codonByIndex(mutation.codonNumber)?.map { it -> it.start()..it.end() } ?: emptyList()
            is GenericRangeMutations -> gene.codingRangeByGenomicCoordinates(mutation.startPosition, mutation.endPosition)
                    ?.map { it -> it.start()..it.end() } ?: emptyList()
        }
    }

    private fun createGeneModel(genericMutations: List<GenericMutation>): Map<String, HmfTranscriptRegion?> {
        val hmfTranscriptRegions = genericMutations.map { geneModel.hmfTranscriptRegionForGenericMutation(it) }
        return hmfTranscriptRegions.filterNotNull().associateBy({ it.gene() }, { it })
    }
}
