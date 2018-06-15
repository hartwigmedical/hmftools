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
        "returns empty for negative codon" {
            forwardGene.codingRanges(-1).size shouldBe 0
            reverseGene.codingRanges(-1).size shouldBe 0
        }

        "returns empty for codon zero" {
            forwardGene.codingRanges(0).size shouldBe 0
            reverseGene.codingRanges(0).size shouldBe 0
        }

        "returns empty for codon five" {
            forwardGene.codingRanges(5).size shouldBe 0
            reverseGene.codingRanges(5).size shouldBe 0
        }

        "returns correct coding ranges for codon in first exon" {
            forwardGene.codingRanges(1) shouldBe listOf(2L..4L)
            reverseGene.codingRanges(1) shouldBe listOf(12L..14L)
            reverseGene.codingRanges(2) shouldBe listOf(9L..11L)
        }

        "returns correct coding ranges for codon in last exon" {
            forwardGene.codingRanges(3) shouldBe listOf(9L..11L)
            forwardGene.codingRanges(4) shouldBe listOf(12L..14L)
            reverseGene.codingRanges(4) shouldBe listOf(2L..4L)
        }

        "returns correct coding ranges for codon overlapping exons" {
            forwardGene.codingRanges(2) shouldBe listOf(5L..5L, 7L..8L)
            reverseGene.codingRanges(3) shouldBe listOf(5L..5L, 7L..8L)
        }

        "returns correct coding ranges for whole gene" {
            forwardGene.codingRanges() shouldBe listOf(2L..5L, 7L..14L)
            reverseGene.codingRanges() shouldBe listOf(2L..5L, 7L..14L)
        }

        "returns correct coding ranges for codon range in exon" {
            forwardGene.codingRanges(3, 4) shouldBe listOf(9L..14L)
            reverseGene.codingRanges(1, 2) shouldBe listOf(9L..14L)
        }

        "returns correct coding ranges for codon range spanning exons" {
            forwardGene.codingRanges(1, 2) shouldBe listOf(2L..5L, 7L..8L)
            forwardGene.codingRanges(1, 3) shouldBe listOf(2L..5L, 7L..11L)
            reverseGene.codingRanges(1, 3) shouldBe listOf(5L..5L, 7L..14L)
            reverseGene.codingRanges(2, 4) shouldBe listOf(2L..5L, 7L..11L)
        }

        "returns empty for invalid codon ranges" {
            forwardGene.codingRanges(1, 9).size shouldBe 0
            forwardGene.codingRanges(0, 3).size shouldBe 0
            reverseGene.codingRanges(1, 9).size shouldBe 0
            reverseGene.codingRanges(0, 3).size shouldBe 0
        }

        "returns correct last codon number" {
            forwardGene.lastCodon shouldBe 4
            reverseGene.lastCodon shouldBe 4
        }
    }
}
