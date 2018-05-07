package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.ActionableFusionOutput
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene

private const val MIN_PARTNER_COUNT = 2

//MIVO: fusions that appear more than MIN_PARTNER_COUNT on the fiveGene side of known fusion pairs
fun inferredPromiscuousFiveGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.knownFusionPairs }
            .distinct()
            .groupBy { it.fiveGene }
            .filter { it.value.size > MIN_PARTNER_COUNT }
            .map { PromiscuousGene(it.key) }
            .toList()
}

//MIVO: fusions that appear more than MIN_PARTNER_COUNT on the threeGene side of known fusion pairs
fun inferredPromiscuousThreeGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.knownFusionPairs }
            .distinct()
            .groupBy { it.threeGene }
            .filter { it.value.size > MIN_PARTNER_COUNT }
            .map { PromiscuousGene(it.key) }
            .toList()
}

//MIVO: genes that appear as promiscuous in knowledgebases and also appear at least once (as fiveGene) with a specific partner in known fusion pairs
fun externalPromiscuousFiveGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.promiscuousGenes }
            .distinct()
            .filter { knowledgebases.flatMap { it.knownFusionPairs }.distinct().map { it.fiveGene }.toSet().contains(it.gene) }
}

//MIVO: genes that appear as promiscuous in knowledgebases and also appear at least once (as threeGene) with a specific partner in known fusion pairs
fun externalPromiscuousThreeGenes(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return knowledgebases.flatMap { it.promiscuousGenes }
            .distinct()
            .filter { knowledgebases.flatMap { it.knownFusionPairs }.distinct().map { it.threeGene }.toSet().contains(it.gene) }
}

fun knownFusionPairs(knowledgebases: Collection<Knowledgebase>): List<FusionPair> {
    return knowledgebases.flatMap { it.knownFusionPairs }.distinct()
}

fun knownPromiscuousFive(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return (inferredPromiscuousFiveGenes(knowledgebases) + externalPromiscuousFiveGenes(knowledgebases)).distinct()
}

fun knownPromiscuousThree(knowledgebases: Collection<Knowledgebase>): List<PromiscuousGene> {
    return (inferredPromiscuousThreeGenes(knowledgebases) + externalPromiscuousThreeGenes(knowledgebases)).distinct()
}

fun actionableFusionPairs(knowledgebases: Collection<Knowledgebase>): List<ActionableFusionOutput> {
    return knowledgebases.flatMap { it.actionableFusions }.filter { it.fusion is FusionPair }.distinct()
}

fun actionablePromiscuousFive(knowledgebases: Collection<Knowledgebase>): List<ActionableFusionOutput> {
    val knownPromiscuousFive = knownPromiscuousFive(knowledgebases).toSet()
    return actionablePromiscuousGenes(knowledgebases).filter { knownPromiscuousFive.contains(it.fusion) }
}

fun actionablePromiscuousThree(knowledgebases: Collection<Knowledgebase>): List<ActionableFusionOutput> {
    val knownPromiscuousThree = knownPromiscuousThree(knowledgebases).toSet()
    return actionablePromiscuousGenes(knowledgebases).filter { knownPromiscuousThree.contains(it.fusion) }
}

private fun actionablePromiscuousGenes(knowledgebases: Collection<Knowledgebase>): List<ActionableFusionOutput> {
    return knowledgebases.flatMap { it.actionableFusions }.filter { it.fusion is PromiscuousGene }
}
