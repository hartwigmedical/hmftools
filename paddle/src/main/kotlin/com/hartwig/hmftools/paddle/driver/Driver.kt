package com.hartwig.hmftools.paddle.driver

import com.hartwig.hmftools.paddle.dnds.DndsCv


data class GeneVariants(val known: Int, val unknown: Double) {
    val total = known + unknown
}

data class DriverLikelihood(private val dndsRate: DndsCv, private val geneVariants: GeneVariants) {
    val expectedDrivers = dndsRate.expectedDrivers(geneVariants.total)
    val driverLikelihood = ((expectedDrivers - geneVariants.known) / geneVariants.unknown).coerceIn(0.0, 1.0)
    val geneUnknownDrivers = geneVariants.unknown * driverLikelihood
    val genePassengers = geneVariants.unknown - geneUnknownDrivers

}
