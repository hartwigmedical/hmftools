package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceMatch


class FragmentAlleles(val aminoAcidFragment: AminoAcidFragment, val full: Collection<HlaAllele>, val partial: Collection<HlaAllele>, val wild: Collection<HlaAllele>) {


    companion object {
        fun create(aminoAcidFragments: List<AminoAcidFragment>, hetLoci: Collection<Int>, sequences: Collection<HlaSequence>, nucleotideLoci: Collection<Int>, nucleotideSequences: Collection<HlaSequence>): List<FragmentAlleles> {
            return aminoAcidFragments.map { create(it, hetLoci, sequences, nucleotideLoci, nucleotideSequences) }.filter { it.full.isNotEmpty() || it.partial.isNotEmpty() }
        }

        private fun create(
                aminoAcidFragment: AminoAcidFragment,
                aminoAcidLoci: Collection<Int>, aminoAcidSequences: Collection<HlaSequence>,
                nucleotideLoci: Collection<Int>, nucleotideSequences: Collection<HlaSequence>): FragmentAlleles {

            val fragmentNucleotideLoci = (aminoAcidFragment.nucleotideIndices() intersect nucleotideLoci).sorted().toIntArray()
            val fragmentNucleotides = fragmentNucleotideLoci.map { aminoAcidFragment.nucleotide(it) }.toCharArray()
            val matchingNucleotideSequences = nucleotideSequences
                    .map { Pair(it.allele, it.match(fragmentNucleotideLoci, fragmentNucleotides)) }
                    .filter { it.second != HlaSequenceMatch.NONE }
                    .filter { aminoAcidFragment.genes.contains("HLA-${it.first.gene}") }
                    .map { Pair(it.first.specificProtein(), it.second) }
                    .distinct()

            val fullNucleotideMatch = matchingNucleotideSequences.filter { it.second == HlaSequenceMatch.FULL }.map { it.first }.toSet()
            val partialNucleotideMatch = matchingNucleotideSequences
                    .filter { it.second == HlaSequenceMatch.PARTIAL  || it.second == HlaSequenceMatch.WILD}
                    .map { it.first }
                    .toSet() subtract fullNucleotideMatch

            val fragmentAminoAcidLoci = (aminoAcidFragment.aminoAcidIndices() intersect aminoAcidLoci).sorted().toIntArray()
            val fragmentAminoAcids = fragmentAminoAcidLoci.map { aminoAcidFragment.aminoAcid(it) }.toCharArray()
            val matchingAminoAcidSequences = aminoAcidSequences
                    .map { Pair(it.allele, it.match(fragmentAminoAcidLoci, fragmentAminoAcids)) }
                    .filter { it.second != HlaSequenceMatch.NONE }
                    .filter { aminoAcidFragment.genes.contains("HLA-${it.first.gene}") }

            val fullAminoAcidMatch = matchingAminoAcidSequences.filter { it.second == HlaSequenceMatch.FULL }.map { it.first }.toSet()
            val partialAminoAcidMatch = matchingAminoAcidSequences.filter { it.second == HlaSequenceMatch.PARTIAL }.map { it.first }.toSet()
            val wildAminoAcidMatch = matchingAminoAcidSequences.filter { it.second == HlaSequenceMatch.WILD }.map { it.first }.toSet()


            if (fullNucleotideMatch.isEmpty() && partialNucleotideMatch.isEmpty()) {
                return FragmentAlleles(aminoAcidFragment, fullAminoAcidMatch, partialAminoAcidMatch, wildAminoAcidMatch)
            }

            val consistentFull = fullAminoAcidMatch.filter { it.specificProtein() in fullNucleotideMatch }
            val remainingFull = fullAminoAcidMatch.filter { it.specificProtein() !in fullNucleotideMatch }

            return FragmentAlleles(aminoAcidFragment, consistentFull, remainingFull union partialAminoAcidMatch, wildAminoAcidMatch)

        }
    }


}