package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.knowledgebaseimporter.Knowledgebase
import com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology.DiseaseOntology
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseSource
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordAnalyzer
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.output.*
import org.apache.logging.log4j.LogManager


class OncoKb(annotatedVariantsLocation: String, actionableVariantsLocation: String, diseaseOntology: DiseaseOntology,
             private val recordAnalyzer: RecordAnalyzer, treatmentTypeMap: Map<String, String>) : Knowledgebase,
        KnowledgebaseSource<OncoKnownRecord, ActionableRecord> {
    private val logger = LogManager.getLogger("OncoKb")

    override val source = "oncoKb"
    override val knownVariants by lazy {
        logger.info("Extracting oncoKb known variants.")
        recordAnalyzer.knownVariants(listOf(this)).distinct()
    }

    override val knownFusionPairs by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<FusionPair>().distinct() }
    override val promiscuousGenes by lazy { knownKbRecords.flatMap { it.events }.filterIsInstance<PromiscuousGene>().distinct() }

    override val actionableVariants by lazy { actionableKbItems.filterIsInstance<ActionableVariantOutput>() }
    override val actionableCNVs by lazy { actionableKbItems.filterIsInstance<ActionableCNVOutput>() }
    override val actionableFusionPairs by lazy { actionableKbItems.filterIsInstance<ActionableFusionPairOutput>() }
    override val actionablePromiscuousGenes by lazy { actionableKbItems.filterIsInstance<ActionablePromiscuousGeneOutput>() }
    override val actionableRanges by lazy { actionableKbItems.filterIsInstance<ActionableGenomicRangeOutput>() }
    override val cancerTypes by lazy {
        actionableKbRecords.flatMap { it.actionability }.map { it.cancerType }
                .associateBy({ it }, { diseaseOntology.findDoids(it) })
    }
    override val knownKbRecords by lazy {
        logger.info("Reading oncoKb known records.")
        CsvReader.readTSVByName<OncoKnownInput>(annotatedVariantsLocation).mapNotNull { it.corrected() }.map { OncoKnownRecord(it) }
    }
    override val actionableKbRecords by lazy {
        logger.info("Reading oncoKb actionable records.")
        CsvReader.readTSVByName<OncoActionableInput>(actionableVariantsLocation, nullString = "null").mapNotNull { it.corrected() }
                .map { OncoActionableRecord(it, treatmentTypeMap) }
    }
    private val actionableKbItems by lazy { recordAnalyzer.actionableItems(listOf(this)).distinct() }
}
