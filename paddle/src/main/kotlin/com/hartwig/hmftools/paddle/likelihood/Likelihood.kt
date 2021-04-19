package com.hartwig.hmftools.paddle.likelihood

import com.google.common.collect.Lists
import com.hartwig.hmftools.paddle.Gene
import com.hartwig.hmftools.paddle.cohort.CohortLoad
import com.hartwig.hmftools.paddle.dnds.DndsCv
import com.hartwig.hmftools.paddle.dnds.DndsCvGene
import com.hartwig.hmftools.paddle.mutation.Mutations
import com.hartwig.hmftools.paddle.mutation.MutationsGene
import java.io.File
import java.nio.file.Files


data class Likelihood(private val cohortSize: Int, private val tumorMutationalLoad: Int, private val dndsRate: DndsCv, private val mutations: Mutations) {
    val expectedDrivers = dndsRate.expectedDrivers(mutations.total)

    val vusDrivers = (expectedDrivers - mutations.known).coerceAtLeast(0.0)
    val vusDriversPerSample = vusDrivers / cohortSize

    val passengers = mutations.unknown - vusDrivers
    val passengersPerMutation = passengers / tumorMutationalLoad

    override fun toString(): String {
        return "$vusDriversPerSample\t$passengersPerMutation"
    }
}

data class LikelihoodGene(
        val gene: Gene,
        val synonymous: Int,
        val redundant: Int,
        val missense: Likelihood,
        val nonsense: Likelihood,
        val splice: Likelihood,
        val indel: Likelihood) {

    companion object {

        operator fun invoke(load: CohortLoad, dndsMap: Map<Gene, DndsCvGene>, mutationsMap: Map<Gene, MutationsGene>): Map<Gene, LikelihoodGene> {
            val result = mutableMapOf<Gene, LikelihoodGene>()

            for (entry in dndsMap.entries) {
                val dnds = entry.value
                mutationsMap[entry.key]?.let { result.put(entry.key, invoke(load, dnds, it)) }
            }

            return result
        }

        operator fun invoke(load: CohortLoad, dnds: DndsCvGene, mutations: MutationsGene): LikelihoodGene {
            if (mutations.gene != dnds.gene) {
                throw IllegalArgumentException("Unable to combine results from difference genes: ${mutations.gene} & ${dnds.gene}")
            }

            val missense = Likelihood(load.cohortSize, load.snv, dnds.missense, mutations.missense)
            val nonsense = Likelihood(load.cohortSize, load.snv, dnds.nonsense, mutations.nonsense)
            val splice = Likelihood(load.cohortSize, load.snv, dnds.splice, mutations.splice)
            val indel = Likelihood(load.cohortSize, load.indel, dnds.indel, mutations.frameshift + mutations.inframe)

            return LikelihoodGene(mutations.gene, mutations.synonymous, mutations.redundant, missense, nonsense, splice, indel)
        }

        fun writeFile(missenseOnly: Boolean, filename: String, likelihoods: Collection<LikelihoodGene>) {
            Files.write(File(filename).toPath(), toLines(missenseOnly, likelihoods))
        }

        private fun toLines(missenseOnly: Boolean, likelihoods: Collection<LikelihoodGene>): List<String> {
            val lines: MutableList<String> = Lists.newArrayList()
            lines.add(headerString(missenseOnly))
            likelihoods.sortedBy { x -> x.gene }.map { it.toString(missenseOnly) }.forEach { e: String -> lines.add(e) }
            return lines
        }

        private fun headerString(missenseOnly: Boolean): String {
            fun headerString(prefix: String): String {
                return "${prefix}VusDriversPerSample\t${prefix}PassengersPerMutation"
            }

            val missenseHeader = "gene\t${headerString("missense")}"
            return if (missenseOnly) missenseHeader else "$missenseHeader\t${headerString("nonsense")}\t${headerString("splice")}\t${headerString("indel")}"
        }
    }

    private fun toString(missenseOnly: Boolean): String {
        return if (missenseOnly) "$gene\t$missense" else "$gene\t$missense\t$nonsense\t$splice\t$indel"
    }
}
