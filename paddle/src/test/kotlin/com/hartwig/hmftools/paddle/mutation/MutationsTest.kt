package com.hartwig.hmftools.paddle.mutation

import com.hartwig.hmftools.paddle.Impact
import com.hartwig.hmftools.paddle.dnds.DndsMutation
import com.hartwig.hmftools.paddle.dnds.DndsMutationTest
import com.hartwig.hmftools.paddle.dnds.DndsMutationTest.Companion.GENE
import org.junit.Assert

object MutationsTest {

    fun mutation(known: Boolean, impact: Impact): DndsMutation {
        return DndsMutationTest.dndsMutation(GENE, known, false, 0, impact)
    }

    fun synonymousMutation(): DndsMutation {
        return DndsMutationTest.dndsMutation(GENE, false, false, 0, Impact.SYNONYMOUS)
    }

    internal fun assertRedundant(expected: Int, victim: MutationsGene) {
        Assert.assertEquals(expected, victim.redundant)
    }

    internal fun assertEmpty(victim: Mutations) {
        Assert.assertEquals(0, victim.known)
        Assert.assertEquals(0, victim.unknown)
    }

    internal fun assertMutations(expectedKnown: Int, expectedUnknown: Int, victim: Mutations) {
        Assert.assertEquals(expectedKnown, victim.known)
        Assert.assertEquals(expectedUnknown, victim.unknown)
        Assert.assertEquals(expectedKnown + expectedUnknown, victim.total)
    }
}