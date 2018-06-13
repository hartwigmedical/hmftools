package com.hartwig.hmftools.knowledgebaseimporter.dao

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class GeneTest : StringSpec() {
    //MIVO:     both forward and reverse cases use the following scenario (* = UTR):
    //exons:        |--------------------|   |------------------------------------------------|
    //positions:    --1*---2---3---4---5---x---7---8---9---10---11---12---13---14---15*---16*--
    //codons:            |___________|_______________|_____________|______________|

    private val forwardStrandExons = listOf(Exon(1, 5, -1, 1), Exon(7, 16, 1, -1))
    private val forwardGene = Gene(forwardStrandExons, 2, 8)

    private val reverseStrandExons = listOf(Exon(-16, -7, -1, 2), Exon(-5, -1, 2, -1))
    private val reverseGene = Gene(reverseStrandExons, 3, 4)

    init {
        "returns null for negative codon" {
            forwardGene.codonPositions(-1) shouldBe null
            reverseGene.codonPositions(-1) shouldBe null
        }

        "returns null for codon zero" {
            forwardGene.codonPositions(0) shouldBe null
            reverseGene.codonPositions(0) shouldBe null
        }

        "returns null for codon five" {
            forwardGene.codonPositions(5) shouldBe null
            reverseGene.codonPositions(5) shouldBe null
        }

        "returns correct positions for codon in first exon" {
            forwardGene.codonPositions(1) shouldBe Triple(2L, 3L, 4L)
            reverseGene.codonPositions(1) shouldBe Triple(12L, 13L, 14L)
            reverseGene.codonPositions(2) shouldBe Triple(9L, 10L, 11L)
        }

        "returns correct positions for codon in last exon" {
            forwardGene.codonPositions(3) shouldBe Triple(9L, 10L, 11L)
            forwardGene.codonPositions(4) shouldBe Triple(12L, 13L, 14L)
            reverseGene.codonPositions(4) shouldBe Triple(2L, 3L, 4L)
        }

        "returns correct positions for codon overlapping exons" {
            forwardGene.codonPositions(2) shouldBe Triple(5L, 7L, 8L)
            reverseGene.codonPositions(3) shouldBe Triple(5L, 7L, 8L)
        }
    }
}
