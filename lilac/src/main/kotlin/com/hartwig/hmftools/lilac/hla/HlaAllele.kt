package com.hartwig.hmftools.lilac.hla

import org.apache.logging.log4j.util.Strings

data class HlaAllele(val gene: String, val alleleGroup: String, val protein: String, val remainder: String) : Comparable<HlaAllele> {

    companion object {
        operator fun invoke(contig: String): HlaAllele {
            val starIndex = contig.indexOf("*")
            val gene = contig.substring(0, starIndex)
            val contigRemainder = contig.substring(starIndex + 1)
            val contigSplit = contigRemainder.split(":")
            val alleleGroup = contigSplit[0]
            val protein = if (contigSplit.size == 1) Strings.EMPTY else contigSplit[1]
            val finalRemainder = if (contigSplit.size == 1) Strings.EMPTY else  contigRemainder.substring(alleleGroup.length + protein.length + 1)
            return HlaAllele(gene, alleleGroup, protein, finalRemainder)
        }
    }

    override fun toString(): String {
        return "${fourDigitName()}$remainder"
    }

    fun fourDigitName(): String {
        return "$gene*$alleleGroup:$protein"
    }

    fun specificProtein(): HlaAllele {
        return HlaAllele(gene, alleleGroup, protein, Strings.EMPTY)
    }

    fun alleleGroup(): HlaAllele {
        return HlaAllele(gene, alleleGroup, Strings.EMPTY, Strings.EMPTY)
    }

    override fun compareTo(other: HlaAllele): Int {
        val geneCompare = gene.compareTo(other.gene)
        if (geneCompare != 0) {
            return geneCompare
        }

        val groupCompare = alleleGroup.compareTo(other.alleleGroup)
        if (groupCompare != 0) {
            return groupCompare
        }

        val proteinCompare = protein.compareTo(other.protein)
        if (proteinCompare != 0) {
            return proteinCompare
        }

        return remainder.compareTo(other.remainder)
    }
}