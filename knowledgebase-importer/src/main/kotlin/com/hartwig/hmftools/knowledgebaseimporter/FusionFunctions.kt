package com.hartwig.hmftools.knowledgebaseimporter

private const val GENE_PATTERN = "[A-Za-z0-9-]"
private const val GENE_GROUP = "($GENE_PATTERN+)"

fun extractFusion(gene: String, fusionString: String, separators: List<String>): Fusion {
    return separators.map { extractFusion(gene, fusionString, it) }
            .sortedBy {
                when (it) {
                    is FusionPair      -> 0
                    is PromiscuousGene -> 1
                }
            }
            .first()
}

fun extractFusion(gene: String, fusionString: String, separator: String): Fusion {
    val fiveGene = fiveGene(gene, fusionString, separator)
    val threeGene = threeGene(gene, fusionString, separator)
    return if (fiveGene == null || threeGene == null) {
        PromiscuousGene(gene)
    } else {
        FusionPair(fiveGene, threeGene)
    }
}

fun fiveGene(gene: String, fusion: String, separator: String): String? {
    return if (isFiveGene(gene, fusion, separator)) gene else extractFiveGene(gene, fusion, separator)
}

fun threeGene(gene: String, fusion: String, separator: String): String? {
    return if (isThreeGene(gene, fusion, separator)) gene else extractThreeGene(gene, fusion, separator)
}

fun extractFiveGene(gene: String, fusion: String, separator: String): String? {
    return extractGene(fusion, "$GENE_GROUP${Regex.escape(separator)}${geneStartLetters(gene)}")
}

fun extractThreeGene(gene: String, fusion: String, separator: String): String? {
    return extractGene(fusion, "${geneStartLetters(gene)}$GENE_PATTERN*${Regex.escape(separator)}$GENE_GROUP")
}

private fun extractGene(fusion: String, patternString: String): String? {
    val matchResult = patternString.toRegex().find(fusion)
    return matchResult?.groupValues?.get(1)
}

fun isFiveGene(gene: String, fusion: String, separator: String): Boolean {
    return !isThreeGene(gene, fusion, separator) && fusion.contains(geneStartLetters(gene))
}

fun isThreeGene(gene: String, fusion: String, separator: String): Boolean {
    return fusion.contains("$separator${geneStartLetters(gene)}")
}

private fun geneStartLetters(gene: String) = gene.substring(0, 3)

fun flipFusion(fusion: Fusion, pairsToFlip: Set<FusionPair>): Fusion {
    return if (fusion is FusionPair && pairsToFlip.contains(fusion)) {
        FusionPair(fusion.threeGene, fusion.fiveGene)
    } else {
        fusion
    }
}
