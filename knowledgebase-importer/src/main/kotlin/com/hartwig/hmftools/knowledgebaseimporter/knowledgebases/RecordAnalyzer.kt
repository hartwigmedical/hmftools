package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import com.hartwig.hmftools.knowledgebaseimporter.dao.EnsemblGeneDAO
import com.hartwig.hmftools.knowledgebaseimporter.dao.Gene
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.GenomicRangeEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.KnownVariantOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.SomaticVariantEvent
import com.hartwig.hmftools.knowledgebaseimporter.transvar.*
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.logging.log4j.LogManager

class RecordAnalyzer(transvarLocation: String, private val reference: IndexedFastaSequenceFile, private val geneDAO: EnsemblGeneDAO) {
    companion object {
        private val logger = LogManager.getLogger("RecordAnalyzer")
        private val blacklistedDrugs = setOf("chemotherapy", "aspirin", "steroids")
    }

    private val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    private val proteinAnalyzer = TransvarProteinAnalyzer(transvarLocation)

    fun knownVariants(knowledgebases: List<KnowledgebaseSource<*, *>>): List<KnownVariantOutput> {
        val records = knowledgebases.flatMap { it.knownKbRecords }
        val knownSomaticVariants = extractSomaticVariants(records)
                .filter { it.second is SomaticVariantEvent }.map { Pair(it.first, it.second as SomaticVariantEvent) }
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

    private fun <R : KnowledgebaseRecord> extractSomaticVariants(records: List<R>): List<Pair<R, ActionableEvent>> {
        return extractGdnaVariants(records) + extractTransvarVariants(records) + extractVariants(records)
    }

    private fun <R : KnowledgebaseRecord> extractGdnaVariants(records: List<R>): List<Pair<R, ActionableEvent>> {
        val recordEventPairs = records.collectEvents<GDnaVariant, R>()
        return recordEventPairs.flatMap { (record, gdnaVariant) ->
            listOfNotNull(extractVariant(record.gene, record.transcript, gdnaVariant.gDnaImpact, reference))
                    .map { Pair(record, it) }
        }
    }

    private fun <R : KnowledgebaseRecord> extractTransvarVariants(records: List<R>): List<Pair<R, ActionableEvent>> {
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
            : List<Pair<R, ActionableEvent>> where T : HgvsAnnotation, R : KnowledgebaseRecord {
        val recordEventPairs = records.collectEvents<T, R>()
        val transvarOutput = analyzer.analyze(recordEventPairs.map { it.second })
        return recordEventPairs.map { it.first }.zip(transvarOutput).flatMap { (record, output) ->
            extractVariants(record.gene, record.transcript, output, reference).map { Pair(record, it) }
        }
    }

    private fun <R : ActionableRecord> analyzeGenomicRanges(records: List<R>): List<Pair<R, GenomicRangeEvent>> {
        val genericMutationRecords = records.collectEvents<GenericMutation, R>()
        val geneModel = createGeneModel(genericMutationRecords.map { it.second })
        return genericMutationRecords.flatMap { (record, mutation) ->
            val gene = geneModel[mutation.transcript] ?: geneModel[mutation.gene]
            mutationCodingRange(mutation, gene).map {
                Pair(record, GenomicRangeEvent(mutation.gene, mutation.transcript.orEmpty(), gene!!.chromosome, it.start.toString(),
                                               it.endInclusive.toString(), gene.transcript))
            }
        }
    }

    private fun mutationCodingRange(mutation: GenericMutation, gene: Gene?): List<ClosedRange<Long>> {
        if (gene == null) {
            logger.warn("Gene model for gene ${mutation.gene}, transcript: ${mutation.transcript} is null")
            return emptyList()
        }
        return when (mutation) {
            is GeneMutations         -> listOf(gene.range())
            is ExonMutations         -> gene.exonCodingRanges(mutation.exonNumber)
            is CodonRangeMutations   -> gene.codonCodingRanges(mutation.startCodon, mutation.endCodon)
            is CodonMutations        -> gene.codonCodingRanges(mutation.codonNumber)
            is GenericRangeMutations -> gene.codingRangesBetween(mutation.startPosition, mutation.endPosition)
        }
    }

    private fun createGeneModel(genericMutations: List<GenericMutation>): Map<String, Gene?> {
        val transcriptSet = genericMutations.mapNotNull { it.transcript }.toSet()
        val geneSet = genericMutations.filter { it.transcript == null }.map { it.gene }.toSet()
        val transcriptGeneModel = transcriptSet.associateBy({ it }, { geneDAO.transcriptGeneModel(it) })
        val canonicalGeneModel = geneSet.associateBy({ it }, { geneDAO.canonicalGeneModel(it) })
        return transcriptGeneModel + canonicalGeneModel
    }
}
