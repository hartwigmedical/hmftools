package com.hartwig.hmftools.paddle.likelihood

import com.hartwig.hmftools.paddle.dnds.DndsCv
import com.hartwig.hmftools.paddle.mutation.Mutations
import org.junit.Assert.assertEquals
import org.junit.Test

class LikelihoodTest {

    private val epsilon = 1e-5

    @Test
    fun testNoVariants() {
        val cohortSize = 3524
        val tumorMutationalLoad = 89924605

        val dndsAR = DndsCv(0, 0.9687)
        val mutations = Mutations(0, 0)
        val victim = Likelihood(cohortSize, tumorMutationalLoad, dndsAR, mutations)
        assertEquals(0.0, victim.vusDriversPerSample, epsilon)
    }

    @Test
    fun testLikelihood() {
        val cohortSize = 3524
        val tumorMutationalLoad = 89924605

        val dndsAR = DndsCv(115, 3.4091126)
        val mutations = Mutations(43, 60)
        val victim = Likelihood(cohortSize, tumorMutationalLoad, dndsAR, mutations)

        assertEquals(72.78686, victim.expectedDrivers, epsilon)

        assertEquals(29.78686, victim.vusDrivers, epsilon)
        assertEquals(30.21314, victim.passengers, epsilon)
        assertEquals(mutations.unknown.toDouble(), victim.vusDrivers + victim.passengers, epsilon)

        assertEquals(8.45257, victim.vusDriversPerSample * 1e3, epsilon)
        assertEquals(3.35983, victim.passengersPerMutation * 1e7, epsilon)
    }
}