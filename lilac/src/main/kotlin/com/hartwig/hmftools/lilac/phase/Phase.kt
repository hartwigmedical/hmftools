package com.hartwig.hmftools.lilac.phase

data class PhasedSequence2(val parent: PhasedSequence?, val index: Int, val contributors: Set<String>, val sequence: List<Pair<Char, Int>>) {

    companion object {

//        fun initialise(index: Int, realigned: List<HlaAlign>) {
//            val candidates = realigned
//                    .filter { it.containsAminoAcid(index) }
//                    .groupBy {  }
//                    .map { Pair(it, it.aminoAcidAt(index, 30)) }
//
//
//            val aminoAcidCounts = candidates.map { it.second }.gr
//
//        }

        fun createNew(index: Int, entry: Char, contributors: Set<String>): PhasedSequence2 {
            return PhasedSequence2(null, index, contributors, listOf(Pair(entry, contributors.size)))
        }


    }

}