package com.hartwig.hmftools.lilackt.hla

import org.apache.logging.log4j.util.Strings

data class HlaAllele(val gene: String, val alleleGroup: String, val protein: String, val synonymous: String, val synonymousNonCoding: String) : Comparable<HlaAllele> {

    companion object {
        operator fun invoke(contig: String): HlaAllele {
            val starIndex = contig.indexOf("*")
            val gene = contig.substring(0, starIndex)
            val contigRemainder = contig.substring(starIndex + 1)
            val contigSplit = contigRemainder.split(":")
            val alleleGroup = contigSplit[0]
            val protein = if (contigSplit.size < 2) Strings.EMPTY else contigSplit[1]
            val synonymousCoding = if (contigSplit.size < 3) Strings.EMPTY else contigSplit[2]
            val synonymousNonCoding = if (contigSplit.size < 4) Strings.EMPTY else contigSplit[3]
            return HlaAllele(gene, alleleGroup, protein, synonymousCoding, synonymousNonCoding)
        }
    }

    override fun toString(): String {
        if (protein.isEmpty()) {
            return "$gene*$alleleGroup"
        }

        if (synonymous.isEmpty()) {
            return "$gene*$alleleGroup:$protein"
        }

        if (synonymousNonCoding.isEmpty()) {
            return "$gene*$alleleGroup:$protein:$synonymous"
        }

        return "$gene*$alleleGroup:$protein:$synonymous:$synonymousNonCoding"
    }

    fun asSixDigit(): HlaAllele {
        return HlaAllele(gene, alleleGroup, protein, synonymous, Strings.EMPTY)
    }

    fun asFourDigit(): HlaAllele {
        return HlaAllele(gene, alleleGroup, protein, Strings.EMPTY, Strings.EMPTY)
    }

    fun asAlleleGroup(): HlaAllele {
        return HlaAllele(gene, alleleGroup, Strings.EMPTY, Strings.EMPTY, Strings.EMPTY)
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

        val synonymousCodingCompare = synonymous.compareTo(other.synonymous)
        if (synonymousCodingCompare != 0) {
            return synonymousCodingCompare
        }

        return synonymousNonCoding.compareTo(other.synonymousNonCoding)
    }
}