package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceMatch


class FragmentSequences(val fragment: Fragment, val full: Collection<HlaAllele>, val partial: Collection<HlaAllele>) {


    companion object {
        fun create(fragments: List<Fragment>, hetLoci: Collection<Int>, sequences: Collection<HlaSequence>): List<FragmentSequences> {
            return fragments.map { create(it, hetLoci, sequences) }.filter { it.full.isNotEmpty() || it.partial.isNotEmpty()}
        }

        fun create(fragment: Fragment, hetLoci: Collection<Int>, sequences: Collection<HlaSequence>): FragmentSequences {
            val loci = (fragment.aminoAcidIndices() intersect hetLoci).sorted().toIntArray()
            val aminoAcids = loci.map { fragment.aminoAcid(it) }.toCharArray()
            val matchingSequences = sequences
                    .map { Pair(it, it.match(loci, aminoAcids)) }
                    .filter { it.second != HlaSequenceMatch.NONE }


            return FragmentSequences(fragment,
                    matchingSequences.filter { it.second == HlaSequenceMatch.FULL }.map { it.first.allele },
                    matchingSequences.filter { it.second == HlaSequenceMatch.PARTIAL }.map { it.first.allele })

        }
    }


}