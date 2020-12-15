package com.hartwig.hmftools.lilac.hla

data class HlaAllele(val gene: String, val alleleGroup: String, val protein: String, val remainder: String) {

    companion object {
        operator fun invoke(contig: String): HlaAllele {
            val starIndex = contig.indexOf("*")
            val gene = contig.substring(0, starIndex)
            val contigRemainder = contig.substring(starIndex + 1)
            val contigSplit = contigRemainder.split(":")
            val alleleGroup = contigSplit[0]
            val protein = contigSplit[1]
            val finalRemainder = contigRemainder.substring(alleleGroup.length + protein.length + 1)
            return HlaAllele(gene, alleleGroup, protein, finalRemainder);
        }
    }

    override fun toString(): String {
        return "${fourDigitName()}$remainder"
    }

    fun fourDigitName(): String {
        return "$gene*$alleleGroup:$protein"
    }
}