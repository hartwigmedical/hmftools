package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableFusionPairOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.ActionablePromiscuousGeneOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene

private const val MIN_PARTNER_COUNT = 2

// Fusions that appear more than MIN_PARTNER_COUNT on the fiveGene side of known fusion pairs
fun inferredPromiscuousFiveGenes(knowledgebases: Knowledgebase): List<PromiscuousGene> {
    return knowledgebases.knownFusionPairs
            .distinct()
            .groupBy { it.fiveGene }
            .filter { it.value.size > MIN_PARTNER_COUNT }
            .map { PromiscuousGene(it.key) }
            .toList()
}

// Fusions that appear more than MIN_PARTNER_COUNT on the fiveGene side of known fusion pairs
fun inferredPromiscuousFiveGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.knownFusionPairs }
            .distinct()
            .groupBy { it.fiveGene }
            .filter { it.value.size > MIN_PARTNER_COUNT }
            .map { PromiscuousGene(it.key) }
            .toList()
}

// Fusions that appear more than MIN_PARTNER_COUNT on the threeGene side of known fusion pairs
fun inferredPromiscuousThreeGenes(knowledgebases: Knowledgebase): List<PromiscuousGene> {
    return knowledgebases.knownFusionPairs
            .distinct()
            .groupBy { it.threeGene }
            .filter { it.value.size > MIN_PARTNER_COUNT }
            .map { PromiscuousGene(it.key) }
            .toList()
}

// Fusions that appear more than MIN_PARTNER_COUNT on the threeGene side of known fusion pairs
fun inferredPromiscuousThreeGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.knownFusionPairs }
            .distinct()
            .groupBy { it.threeGene }
            .filter { it.value.size > MIN_PARTNER_COUNT }
            .map { PromiscuousGene(it.key) }
            .toList()
}

// Genes that appear as promiscuous in knowledgebases and also appear at least once (as fiveGene) with a specific partner in known fusion pairs
fun externalPromiscuousFiveGenes(knowledgebases: Knowledgebase): List<PromiscuousGene> {
    return knowledgebases.promiscuousGenes
            .distinct()
            .filter { knowledgebases.knownFusionPairs.distinct().map { it.fiveGene }.toSet().contains(it.gene) }
}

// Genes that appear as promiscuous in knowledgebases and also appear at least once (as fiveGene) with a specific partner in known fusion pairs
fun externalPromiscuousFiveGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.promiscuousGenes }
            .distinct()
            .filter { knowledgebases.flatMap { it.knownFusionPairs }.distinct().map { it.fiveGene }.toSet().contains(it.gene) }
}

// Genes that appear as promiscuous in knowledgebases and also appear at least once (as threeGene) with a specific partner in known fusion pairs
fun externalPromiscuousThreeGenes(knowledgebases: Knowledgebase): List<PromiscuousGene> {
    return knowledgebases.promiscuousGenes
            .distinct()
            .filter { knowledgebases.knownFusionPairs.distinct().map { it.threeGene }.toSet().contains(it.gene) }
}

// Genes that appear as promiscuous in knowledgebases and also appear at least once (as threeGene) with a specific partner in known fusion pairs
fun externalPromiscuousThreeGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.promiscuousGenes }
            .distinct()
            .filter { knowledgebases.flatMap { it.knownFusionPairs }.distinct().map { it.threeGene }.toSet().contains(it.gene) }
}

fun knownFusionPairs(knowledgebases: Knowledgebase): List<FusionPair> {
    return knowledgebases.knownFusionPairs.distinct()
}

fun knownPromiscuousFive(knowledgebases: Knowledgebase): List<PromiscuousGene> {
    return (inferredPromiscuousFiveGenes(knowledgebases) + externalPromiscuousFiveGenes(knowledgebases)).distinct()
}

fun knownPromiscuousFive(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return (inferredPromiscuousFiveGenes(knowledgebases) + externalPromiscuousFiveGenes(knowledgebases)).distinct()
}

fun knownPromiscuousThree(knowledgebases: Knowledgebase): List<PromiscuousGene> {
    return (inferredPromiscuousThreeGenes(knowledgebases) + externalPromiscuousThreeGenes(knowledgebases)).distinct()
}

fun knownPromiscuousThree(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return (inferredPromiscuousThreeGenes(knowledgebases) + externalPromiscuousThreeGenes(knowledgebases)).distinct()
}

fun actionableFusionPairs(knowledgebases: Collection<Knowledgebase>): List<ActionableFusionPairOutput> {
    return knowledgebases.flatMap { it.actionableFusionPairs }.distinct()
}

fun actionablePromiscuousFive(knowledgebases: Collection<Knowledgebase>): List<ActionablePromiscuousGeneOutput> {
    val knownPromiscuousFive = knownPromiscuousFive(knowledgebases).toSet()
    return actionablePromiscuousGenes(knowledgebases).filter { knownPromiscuousFive.contains(it.event) }
}

fun actionablePromiscuousThree(knowledgebases: Collection<Knowledgebase>): List<ActionablePromiscuousGeneOutput> {
    val knownPromiscuousThree = knownPromiscuousThree(knowledgebases).toSet()
    return actionablePromiscuousGenes(knowledgebases).filter { knownPromiscuousThree.contains(it.event) }
}

private fun actionablePromiscuousGenes(knowledgebases: Collection<Knowledgebase>): List<ActionablePromiscuousGeneOutput> {
    return knowledgebases.flatMap { it.actionablePromiscuousGenes }.distinct()
}
