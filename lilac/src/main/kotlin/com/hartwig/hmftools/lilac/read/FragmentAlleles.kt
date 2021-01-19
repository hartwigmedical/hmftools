package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceMatch


class FragmentAlleles(val fragment: Fragment, val full: Collection<HlaAllele>, val partial: Collection<HlaAllele>) {


    companion object {
        fun create(fragments: List<Fragment>, hetLoci: Collection<Int>, sequences: Collection<HlaSequence>, nucleotideLoci: Collection<Int>, nucleotideSequences: Collection<HlaSequence>): List<FragmentAlleles> {
            return fragments.map { create(it, hetLoci, sequences, nucleotideLoci, nucleotideSequences) }.filter { it.full.isNotEmpty() || it.partial.isNotEmpty() }
        }

        private fun create(
                fragment: Fragment,
                aminoAcidLoci: Collection<Int>, aminoAcidSequences: Collection<HlaSequence>,
                nucleotideLoci: Collection<Int>, nucleotideSequences: Collection<HlaSequence>): FragmentAlleles {

            val fragmentNucleotideLoci = (fragment.nucleotideIndices() intersect nucleotideLoci).sorted().toIntArray()
            val fragmentNucleotides = fragmentNucleotideLoci.map { fragment.nucleotide(it) }.toCharArray()
            val matchingNucleotideSequences = nucleotideSequences
                    .map { Pair(it.allele, it.match(fragmentNucleotideLoci, fragmentNucleotides)) }
                    .filter { it.second != HlaSequenceMatch.NONE }
                    .filter { fragment.genes.contains("HLA-${it.first.gene}") }
                    .map { Pair(it.first.specificProtein(), it.second) }
                    .distinct()

            val fullNucleotideMatch = matchingNucleotideSequences.filter { it.second == HlaSequenceMatch.FULL }.map { it.first }.toSet()
            val partialNucleotideMatch = matchingNucleotideSequences.filter { it.second == HlaSequenceMatch.PARTIAL }.map { it.first }.toSet() subtract fullNucleotideMatch


            val fragmentAminoAcidLoci = (fragment.aminoAcidIndices() intersect aminoAcidLoci).sorted().toIntArray()
            val fragmentAminoAcids = fragmentAminoAcidLoci.map { fragment.aminoAcid(it) }.toCharArray()
            val matchingAminoAcidSequences = aminoAcidSequences
                    .map { Pair(it.allele, it.match(fragmentAminoAcidLoci, fragmentAminoAcids)) }
                    .filter { it.second != HlaSequenceMatch.NONE }
                    .filter { fragment.genes.contains("HLA-${it.first.gene}") }

            val fullAminoAcidMatch = matchingAminoAcidSequences.filter { it.second == HlaSequenceMatch.FULL }.map { it.first }.toSet()
            val partialAminoAcidMatch = matchingAminoAcidSequences.filter { it.second == HlaSequenceMatch.PARTIAL }.map { it.first }.toSet()

//            if (fullAminoAcidMatch.size == 1 && fullAminoAcidMatch.contains(HlaAllele("C*07:57"))) {
//                println(fragment.id + " [" + fragmentAminoAcidLoci.joinToString(", ") + "] [" + fragmentAminoAcids.joinToString(", ") + "]")
//            }
//
//            if(fragment.id == "A00624:8:HHKYHDSXX:1:1549:11858:36307 -> 29911268") {
//                println("sdf")
//            }

            if (fullNucleotideMatch.isEmpty() && partialNucleotideMatch.isEmpty()) {
                return FragmentAlleles(fragment, fullAminoAcidMatch, partialAminoAcidMatch)
            }

            val consistentFull = fullAminoAcidMatch.filter { it.specificProtein() in fullNucleotideMatch }
            val remainingFull = fullAminoAcidMatch.filter { it.specificProtein() in partialNucleotideMatch }
            val consistentPartial = partialAminoAcidMatch.filter { it.specificProtein() in fullNucleotideMatch || it.specificProtein() in partialNucleotideMatch }



            if (fragmentNucleotideLoci.contains(1012) || fragmentNucleotideLoci.contains(1013)) {
//                println("Sdf")
            }

            if (!consistentFull.containsAll(fullAminoAcidMatch)) {
//                println("Sdf")
            }
            return FragmentAlleles(fragment, consistentFull, remainingFull union consistentPartial)
//            return FragmentAlleles(fragment, fullAminoAcidMatch, partialAminoAcidMatch)

        }
    }


}