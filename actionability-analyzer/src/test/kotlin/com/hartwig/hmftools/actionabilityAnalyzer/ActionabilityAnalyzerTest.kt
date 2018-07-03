package com.hartwig.hmftools.actionabilityAnalyzer

import com.google.common.io.Resources
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

private const val ON_LABEL = "On-label"

class ActionabilityAnalyzerTest : StringSpec() {
    private val samplesMap = mapOf("CPCT99110022T" to "Skin", "CPCT99110033T" to "Lung")

    private val actionableCNVs = Resources.getResource("actionableCNVs").path
    private val actionableFusionPairs = Resources.getResource("actionableFusionPairs").path
    private val actionablePromiscuousFive = Resources.getResource("actionablePromiscuousFive").path
    private val actionablePromiscuousThree = Resources.getResource("actionablePromiscuousThree").path
    private val actionableVariants = Resources.getResource("actionableVariants").path
    private val actionableRanges = Resources.getResource("actionableRanges").path
    private val cancerTypesMapping = Resources.getResource("knowledgebaseCancerTypes").path
    private val actionabilityAnalyzer = ActionabilityAnalyzer(samplesMap, actionableVariants, actionableFusionPairs,
                                                              actionablePromiscuousFive, actionablePromiscuousThree, actionableCNVs,
                                                              cancerTypesMapping, actionableRanges)

    private val brafSNV = CohortMutation("CPCT99110022T", "7", "140453136", "A", "T", "SNP", "BRAF",
                                         "missense", "missense", "missense", "ENST00000288602", "TRUE")
    private val brafOtherSNV = CohortMutation("CPCT99110033T", "7", "140453136", "A", "T", "SNP", "BRAF",
                                              "missense", "missense", "missense", "ENST00000288602", "TRUE")

    private val ptenSNV = CohortMutation("CPCT99110033T", "10", "89653781", "G", "C", "SNP", "PTEN",
                                         "splice;intron", "splice;intron", "splice;intron", "ENST00000371953", "FALSE")

    init {
        "finds BRAF SNV actionability" {
            val actionability = actionabilityAnalyzer.actionabilityForVariant(brafSNV)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            actionability.size shouldBe 6
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 2
            events.size shouldBe 1
            events.first() shouldBe "${brafSNV.gene} ${brafSNV.chromosome}:${brafSNV.position} ${brafSNV.ref}->${brafSNV.alt}"
            (drugs == setOf("Dabrafenib", "Vemurafenib")) shouldBe true
            (sources == setOf("civic", "cgi")) shouldBe true
        }

        "finding variant does not depend on gene name" {
            val variant = brafSNV.copy(gene = "BRAF2")
            val brafSnvActionability = actionabilityAnalyzer.actionabilityForVariant(variant)
            val actionability = actionabilityAnalyzer.actionabilityForVariant(variant)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            val brafEvents = brafSnvActionability.map { it.event }.toSet()
            actionability.size shouldBe 6
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 2
            events.size shouldBe 1
            events shouldBe brafEvents
            (drugs == setOf("Dabrafenib", "Vemurafenib")) shouldBe true
            (sources == setOf("civic", "cgi")) shouldBe true
        }

        "finds BRAF SNV actionability for different primary tumor" {
            val variant = brafOtherSNV
            val actionability = actionabilityAnalyzer.actionabilityForVariant(variant)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            actionability.size shouldBe 6
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 3
            events.size shouldBe 1
            events.first() shouldBe "${variant.gene} ${variant.chromosome}:${variant.position} ${variant.ref}->${variant.alt}"
            (drugs == setOf("Dabrafenib", "Vemurafenib")) shouldBe true
            (sources == setOf("oncoKb", "cgi")) shouldBe true
        }

        "does not find SNV with different chromosome" {
            actionabilityAnalyzer.actionabilityForVariant(brafSNV.copy(chromosome = "8")).size shouldBe 0
        }

        "does not find SNV with different position" {
            actionabilityAnalyzer.actionabilityForVariant(brafSNV.copy(position = "140453137")).size shouldBe 0
        }

        "does not find SNV with different ref" {
            actionabilityAnalyzer.actionabilityForVariant(brafSNV.copy(ref = "C")).size shouldBe 0
        }

        "does not find SNV with different alt" {
            actionabilityAnalyzer.actionabilityForVariant(brafSNV.copy(alt = "G")).size shouldBe 0
        }

        "finds BRAF SNV range actionability" {
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(brafOtherSNV)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            actionability.size shouldBe 3
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 2
            events.size shouldBe 1
            (drugs == setOf("Vemurafenib", "Dabrafenib")) shouldBe true
            (sources == setOf("civic")) shouldBe true
        }

        "finds BRAF actionability that intersects range" {
            val variant = brafOtherSNV.copy(position = "140453135", ref = "GA", alt = "CT")
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(variant)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            actionability.size shouldBe 3
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 2
            events.size shouldBe 1
            (drugs == setOf("Vemurafenib", "Dabrafenib")) shouldBe true
            (sources == setOf("civic")) shouldBe true
        }

        "finds BRAF actionability that completely overlaps range" {
            val variant = brafOtherSNV.copy(position = "140453135", ref = "GATC", alt = "CTAG")
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(variant)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            actionability.size shouldBe 3
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 2
            events.size shouldBe 1
            (drugs == setOf("Vemurafenib", "Dabrafenib")) shouldBe true
            (sources == setOf("civic")) shouldBe true
        }

        "does not find BRAF actionability before range" {
            val variant = brafOtherSNV.copy(position = "140453135", ref = "G", alt = "C")
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(variant)
            actionability.size shouldBe 0
        }

        "does not find BRAF actionability after range" {
            val variant = brafOtherSNV.copy(position = "140453138", ref = "G", alt = "C")
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(variant)
            actionability.size shouldBe 0
        }

        "does not find BRAF actionability for insert before range" {
            val variant = brafOtherSNV.copy(position = "140453135", ref = "G", alt = "GATC")
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(variant)
            actionability.size shouldBe 0
        }

        "finds BRAF actionability for deletion that intersects range" {
            val variant = brafOtherSNV.copy(position = "140453135", ref = "GA", alt = "G")
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(variant)
            val drugs = drugs(actionability)
            val sources = sources(actionability)
            val events = actionability.map { it.event }.toSet()
            actionability.size shouldBe 3
            actionability.filter { it.treatmentType == ON_LABEL }.size shouldBe 2
            events.size shouldBe 1
            (drugs == setOf("Vemurafenib", "Dabrafenib")) shouldBe true
            (sources == setOf("civic")) shouldBe true
        }

        "does not find PTEN SNV range actionability" {
            val actionability = actionabilityAnalyzer.rangeActionabilityForVariant(ptenSNV)
            actionability.size shouldBe 0
        }
    }

    private fun drugs(actionability: Set<ActionabilityOutput>): Set<String> {
        return actionability.filter { it.treatmentType == ON_LABEL }.map { it.drug }.toSet()
    }

    private fun sources(actionability: Set<ActionabilityOutput>): Set<String> {
        return actionability.filter { it.treatmentType == ON_LABEL }.map { it.source }.toSet()
    }
}
