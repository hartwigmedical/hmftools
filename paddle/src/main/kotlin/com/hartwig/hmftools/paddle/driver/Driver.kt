package com.hartwig.hmftools.paddle.driver

import com.hartwig.hmftools.paddle.Gene
import com.hartwig.hmftools.paddle.cohort.TumorMutationalLoad
import com.hartwig.hmftools.paddle.dnds.DndsCv
import com.hartwig.hmftools.paddle.dnds.DndsCvGene
import com.hartwig.hmftools.paddle.mutation.Mutations
import com.hartwig.hmftools.paddle.mutation.MutationsGene


data class Likelihood(private val cohortSize: Int, private val tumorMutationalLoad: Int, private val dndsRate: DndsCv, private val geneVariants: Mutations) {
    val expectedDrivers = dndsRate.expectedDrivers(geneVariants.total)
    val driverLikelihood = ((expectedDrivers - geneVariants.known) / geneVariants.unknown).coerceIn(0.0, 1.0)
    val geneUnknownDrivers = geneVariants.unknown * driverLikelihood
    val genePassengers = geneVariants.unknown - geneUnknownDrivers
    val pUnknownVariantIsDriver = geneUnknownDrivers / cohortSize
    val expectedPassengerRate = genePassengers / tumorMutationalLoad

    override fun toString(): String {

        //  gene  impact   knownDrivers unknownDrivers expectedDrivers driverLikelihood gene_drivers gene_non_drivers
        //  <chr> <chr>           <dbl>          <dbl>           <dbl>            <dbl>        <dbl>            <dbl>
        //1 AR    Inframe            12             17             0              0              0               17
        //2 AR    Missense           43             60            72.8            0.496         29.8             30.2



        return "(known=${geneVariants.known}, unknown=${geneVariants.unknown}, expectedDrivers=$expectedDrivers)"
//        return "(lle=$driverLikelihood, pDriver=$pUnknownVariantIsDriver, expPassenger=$expectedPassengerRate)"
    }
}

data class LikelihoodGene(val gene: Gene, val missense: Likelihood, val nonsense: Likelihood, val splice: Likelihood, val indel: Likelihood) {
    companion object {

        operator fun invoke(cohortSize: Int, tumorMutationalLoad: TumorMutationalLoad, dndsMap: Map<Gene, DndsCvGene>, mutationsMap: Map<Gene, MutationsGene>): Map<Gene, LikelihoodGene> {
            val result = mutableMapOf<Gene, LikelihoodGene>()

            for (entry in dndsMap.entries) {
                val dnds = entry.value
                mutationsMap[entry.key]?.let { result.put(entry.key, invoke(cohortSize, tumorMutationalLoad, dnds, it)) }
            }

            return result
        }

        operator fun invoke(cohortSize: Int, tumorMutationalLoad: TumorMutationalLoad, dnds: DndsCvGene, mutations: MutationsGene): LikelihoodGene {
            if (mutations.gene != dnds.gene) {
                throw IllegalArgumentException("Unable to combine results from difference genes: ${mutations.gene} & ${dnds.gene}")
            }

            val missense = Likelihood(cohortSize, tumorMutationalLoad.snv, dnds.missense, mutations.missense)
            val nonsense = Likelihood(cohortSize, tumorMutationalLoad.snv, dnds.nonsense, mutations.nonsense)
            val splice = Likelihood(cohortSize, tumorMutationalLoad.snv, dnds.splice, mutations.splice)
            val indel = Likelihood(cohortSize, tumorMutationalLoad.indel, dnds.indel, mutations.frameshift + mutations.inframe)

            return LikelihoodGene(dnds.gene, missense, nonsense, splice, indel)
        }
    }
}
