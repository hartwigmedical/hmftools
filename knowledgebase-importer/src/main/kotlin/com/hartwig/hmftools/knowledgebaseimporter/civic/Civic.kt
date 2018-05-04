package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.ArrayListMultimap
import com.google.common.collect.Multimap
import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.*
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarCdnaAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.transvar.TransvarOutput
import com.hartwig.hmftools.knowledgebaseimporter.transvar.annotations.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.transvar.extractVariants
import com.hartwig.hmftools.knowledgebaseimporter.transvar.somaticVariant
import htsjdk.samtools.reference.IndexedFastaSequenceFile

class Civic(variantsLocation: String, evidenceLocation: String, transvarLocation: String, diseaseOntology: DiseaseOntology,
            private val reference: IndexedFastaSequenceFile) :
        Knowledgebase {
    companion object {
        private const val SOURCE = "civic"
        private val FUSION_SEPARATORS = listOf("-")
        private val FUSIONS_TO_FILTER = setOf(FusionPair("BRAF", "CUL1"))
    }

    private val cdnaAnalyzer = TransvarCdnaAnalyzer(transvarLocation)
    private val records by lazy { preProcessCivicRecords(variantsLocation, evidenceLocation) }
    private val civicVariants by lazy { readCivicVariants() }
    val cancerTypes by lazy { readCancerTypes(diseaseOntology) }

    override val knownVariants: List<KnownVariantOutput> by lazy { knownVariants() }
    override val knownFusionPairs: List<FusionPair> by lazy { actionableFusions.map { it.fusion }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes: List<PromiscuousGene> by lazy { actionableFusions.map { it.fusion }.filterIsInstance<PromiscuousGene>().distinct() }
    override val actionableVariants: List<ActionableVariantOutput> by lazy { actionableVariants() }
    override val actionableCNVs: List<ActionableCNVOutput> by lazy { actionableCNVs() }
    override val actionableFusions: List<ActionableFusionOutput> by lazy { actionableFusions() }


    private fun knownVariants(): List<KnownVariantOutput> {
        return civicVariants.map { (civicRecord, somaticVariant) ->
            KnownVariantOutput(civicRecord.gene, civicRecord.transcript, additionalInfo(civicRecord), somaticVariant)
        }
    }

    private fun actionableVariants(): List<ActionableVariantOutput> {
        return civicVariants.flatMap { (record, somaticVariant) ->
            record.evidence.filter { it.direction == "Supports" }.flatMap { evidence ->
                evidence.drugs.map { drug ->
                    ActionableVariantOutput(record.gene, somaticVariant, actionability(drug, evidence))
                }
            }
        }
    }

    private fun actionableCNVs(): List<ActionableCNVOutput> {
        return records.filter { it.variant == "AMPLIFICATION" || it.variant == "DELETION" || it.variant == "LOH" }.flatMap { record ->
            record.evidence.flatMap { evidence ->
                evidence.drugs.map { drug ->
                    ActionableCNVOutput(record.gene, extractAmpOrDel(record.variant), actionability(drug, evidence))
                }
            }
        }
    }

    private fun actionableFusions(): List<ActionableFusionOutput> {
        val fusionRecords = records.filter { it.variantTypes.contains("fusion") }
        val fusions = fusionRecords.map { extractFusion(it.gene, it.variant.trim(), FUSION_SEPARATORS) }
        return fusionRecords.zip(fusions).filterNot { FUSIONS_TO_FILTER.contains(it.second) }
                .flatMap { (record, fusion) ->
                    record.evidence.flatMap { evidence ->
                        evidence.drugs.map { drug ->
                            ActionableFusionOutput(fusion, actionability(drug, evidence))
                        }
                    }
                }
    }

    private fun readCivicVariants(): List<Pair<CivicRecord, SomaticVariant>> {
        val variantRecords = records.filter { hasVariant(it) }
        val transvarOutput = transvarInferredCdna(variantRecords)
        return variantRecords.zip(transvarOutput)
                .flatMap { (civicRecord, transvarOutput) ->
                    val civicVariant = annotateVariant(civicRecord, reference)
                    val inferredVariants = extractVariants(transvarOutput, reference)
                    (listOfNotNull(civicVariant) + inferredVariants.filterNot { it == civicVariant }).map { Pair(civicRecord, it) }
                }
    }

    private fun transvarInferredCdna(records: List<CivicRecord>): List<TransvarOutput> {
        val cdnaRecords = records.map { record ->
            val hgvsParts = record.hgvs.split(":")
            if (hgvsParts.size > 1) {
                CDnaAnnotation(hgvsParts[0], hgvsParts[1])
            } else {
                CDnaAnnotation("na", "na")
            }
        }
        return cdnaAnalyzer.analyze(cdnaRecords)
    }

    private fun hasVariant(record: CivicRecord): Boolean {
        return (!record.chromosome.isBlank() && !record.start.isBlank() && !record.stop.isBlank() && (!record.ref.isBlank() || !record.alt.isBlank()))
                || (!record.hgvs.isBlank())
    }

    private fun preProcessCivicRecords(variantFileLocation: String, evidenceFileLocation: String): List<CivicRecord> {
        val variantEvidenceMap = readEvidenceMap(evidenceFileLocation)
        return readTSVRecords(variantFileLocation) { CivicRecord(it, variantEvidenceMap) }.map { correctCivicVariants(it) }
    }

    private fun readEvidenceMap(evidenceLocation: String): Multimap<String, CivicEvidence> {
        val evidenceMap = ArrayListMultimap.create<String, CivicEvidence>()
        readTSVRecords(evidenceLocation) { csvRecord -> evidenceMap.put(csvRecord["variant_id"], CivicEvidence(csvRecord)) }
        return evidenceMap
    }

    private fun additionalInfo(civicRecord: CivicRecord): String {
        val highestEvidenceLevel = civicRecord.evidence.map { it.level }.sorted().firstOrNull() ?: "N"
        return (highestEvidenceLevel == "A" || highestEvidenceLevel == "B" || highestEvidenceLevel == "C").toString()
    }

    private fun actionability(drug: String, evidence: CivicEvidence): Actionability {
        return Actionability(SOURCE, evidence.cancerType, drug, evidence.level, evidence.significance, evidence.type)
    }

    private fun extractAmpOrDel(variantType: String): String = if (variantType == "AMPLIFICATION") "Amplification" else "Deletion"

    private fun annotateVariant(civicRecord: CivicRecord, reference: IndexedFastaSequenceFile): SomaticVariant? {
        return if (!civicRecord.chromosome.isEmpty() && !civicRecord.start.isEmpty() && (!civicRecord.ref.isEmpty() || !civicRecord.alt.isEmpty())) {
            val position = civicRecord.start.toLong()
            when {
                civicRecord.ref.isEmpty() -> {
                    val base = reference.getSubsequenceAt(civicRecord.chromosome, position, position).baseString
                    somaticVariant(civicRecord.chromosome, position, base, base + civicRecord.alt)
                }
                civicRecord.alt.isEmpty() -> {
                    val base = reference.getSubsequenceAt(civicRecord.chromosome, position - 1, position - 1).baseString
                    somaticVariant(civicRecord.chromosome, position - 1, base + civicRecord.ref, base)
                }
                else                      -> somaticVariant(civicRecord.chromosome, position, civicRecord.ref, civicRecord.alt)
            }
        } else {
            null
        }
    }

    private fun correctCivicVariants(record: CivicRecord): CivicRecord {
        return when {
            record.variant.contains(Regex("MLL-MLLT3")) && record.gene == "KMT2A" ->
                record.copy(variant = record.variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
            else                                                                  -> record
        }
    }

    private fun readCancerTypes(diseaseOntology: DiseaseOntology): Map<String, Set<Int>> {
        return records.flatMap { it.evidence }.map { Pair(it.cancerType, doidsForEvidence(it, diseaseOntology)) }.toMap()
    }

    private fun doidsForEvidence(evidence: CivicEvidence, diseaseOntology: DiseaseOntology): Set<Int> {
        return diseaseOntology.findDoids(evidence.cancerType) + diseaseOntology.findDoids(evidence.doid.toIntOrNull())
    }
}
