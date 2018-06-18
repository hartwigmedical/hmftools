package com.hartwig.hmftools.knowledgebaseimporter.dao

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class GeneTest : StringSpec() {
    //MIVO:     both forward and reverse cases use the following scenario (* = UTR):
    //exons:        |--------------------|   |------------------------------------------------|
    //positions:    --1*---2---3---4---5---x---7---8---9---10---11---12---13---14---15*---16*--
    //codons:            |___________|_______________|_____________|______________|

    private val forwardStrandExons = listOf(Exon("1", 1, 5, -1, 1),
                                            Exon("1", 7, 16, 1, -1))
    private val forwardGene = Gene(forwardStrandExons, 2, 8)

    private val reverseStrandExons = listOf(Exon("1", -16, -7, -1, 2),
                                            Exon("1", -5, -1, 2, -1))
    private val reverseGene = Gene(reverseStrandExons, 3, 4)

    init {
        "returns empty for negative codon" {
            forwardGene.codonCodingRanges(-1).size shouldBe 0
            reverseGene.codonCodingRanges(-1).size shouldBe 0
        }

        "returns empty for codon zero" {
            forwardGene.codonCodingRanges(0).size shouldBe 0
            reverseGene.codonCodingRanges(0).size shouldBe 0
        }

        "returns empty for codon five" {
            forwardGene.codonCodingRanges(5).size shouldBe 0
            reverseGene.codonCodingRanges(5).size shouldBe 0
        }

        "returns correct coding ranges for codon in first exon" {
            forwardGene.codonCodingRanges(1) shouldBe listOf(2L..4L)
            reverseGene.codonCodingRanges(1) shouldBe listOf(12L..14L)
            reverseGene.codonCodingRanges(2) shouldBe listOf(9L..11L)
        }

        "returns correct coding ranges for codon in last exon" {
            forwardGene.codonCodingRanges(3) shouldBe listOf(9L..11L)
            forwardGene.codonCodingRanges(4) shouldBe listOf(12L..14L)
            reverseGene.codonCodingRanges(4) shouldBe listOf(2L..4L)
        }

        "returns correct coding ranges for codon overlapping exons" {
            forwardGene.codonCodingRanges(2) shouldBe listOf(5L..5L, 7L..8L)
            reverseGene.codonCodingRanges(3) shouldBe listOf(5L..5L, 7L..8L)
        }

        "returns correct coding ranges for whole gene" {
            forwardGene.codingRanges() shouldBe listOf(2L..5L, 7L..14L)
            reverseGene.codingRanges() shouldBe listOf(2L..5L, 7L..14L)
        }

        "returns correct coding ranges for codon range in exon" {
            forwardGene.codonCodingRanges(3, 4) shouldBe listOf(9L..14L)
            reverseGene.codonCodingRanges(1, 2) shouldBe listOf(9L..14L)
        }

        "returns correct coding ranges for codon range spanning exons" {
            forwardGene.codonCodingRanges(1, 2) shouldBe listOf(2L..5L, 7L..8L)
            forwardGene.codonCodingRanges(1, 3) shouldBe listOf(2L..5L, 7L..11L)
            reverseGene.codonCodingRanges(1, 3) shouldBe listOf(5L..5L, 7L..14L)
            reverseGene.codonCodingRanges(2, 4) shouldBe listOf(2L..5L, 7L..11L)
        }

        "returns empty for invalid codon ranges" {
            forwardGene.codonCodingRanges(1, 9).size shouldBe 0
            forwardGene.codonCodingRanges(0, 3).size shouldBe 0
            reverseGene.codonCodingRanges(1, 9).size shouldBe 0
            reverseGene.codonCodingRanges(0, 3).size shouldBe 0
        }

        "returns correct exon coding ranges" {
            forwardGene.exonCodingRanges(1) shouldBe listOf(2L..5L)
            reverseGene.exonCodingRanges(1) shouldBe listOf(7L..14L)
        }

        "returns correct exon coding ranges" {
            forwardGene.exonCodingRanges(1) shouldBe listOf(2L..5L)
            forwardGene.exonCodingRanges(2) shouldBe listOf(7L..14L)
            reverseGene.exonCodingRanges(1) shouldBe listOf(7L..14L)
            reverseGene.exonCodingRanges(2) shouldBe listOf(2L..5L)
        }

        "returns empty for non-existing exons" {
            forwardGene.exonCodingRanges(0).size shouldBe 0
            reverseGene.exonCodingRanges(3).size shouldBe 0
        }

        "returns correct last codon number" {
            forwardGene.lastCodon shouldBe 4
            reverseGene.lastCodon shouldBe 4
        }
    }
}
