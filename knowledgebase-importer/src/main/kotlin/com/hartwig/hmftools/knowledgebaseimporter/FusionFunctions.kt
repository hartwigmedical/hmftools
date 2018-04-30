package com.hartwig.hmftools.knowledgebaseimporter

import org.apache.logging.log4j.LogManager

private val logger = LogManager.getLogger("FusionFunctions")

fun extractFusion(gene: String, fusion: String, separators: List<String>): Pair<String, String>? {
    val fusionPair = separators.mapNotNull { extractFusion(gene, fusion, it) }.firstOrNull()
    if (fusionPair == null) {
        logger.warn("Could not extract fusion pair from $fusion. Searching for gene $gene with separators $separators")
    }
    return fusionPair
}

fun extractFusion(gene: String, fusion: String, separator: String): Pair<String, String>? {
    return if (!isThreeGene(gene, fusion, separator) && !isFiveGene(gene, fusion, separator)) {
        null
    } else {
        Pair(fiveGene(gene, fusion, separator), threeGene(gene, fusion, separator))
    }
}

fun fiveGene(gene: String, fusion: String, separator: String): String {
    return if (isFiveGene(gene, fusion, separator)) gene else fusion.substringBefore("$separator$gene")
}

fun threeGene(gene: String, fusion: String, separator: String): String {
    return if (isThreeGene(gene, fusion, separator)) gene else fusion.substringAfter("$gene$separator")
}

fun isThreeGene(gene: String, fusion: String, separator: String): Boolean {
    return fusion.contains("$separator$gene")
}

fun isFiveGene(gene: String, fusion: String, separator: String): Boolean {
    return fusion.contains("$gene$separator")
}

fun flipGenePair(pair: Pair<String, String>, pairsToFlip: Set<Pair<String, String>>): Pair<String, String> {
    return if (pairsToFlip.contains(pair)) {
        Pair(pair.second, pair.first)
    } else {
        pair
    }
}
