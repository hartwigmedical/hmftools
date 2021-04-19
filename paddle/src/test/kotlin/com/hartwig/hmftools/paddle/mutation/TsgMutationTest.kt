package com.hartwig.hmftools.paddle.mutation

import com.hartwig.hmftools.paddle.Impact
import com.hartwig.hmftools.paddle.dnds.DndsMutationTest
import org.junit.Test

class TsgMutationTest {

    @Test
    fun testSingleUnknownHit() {
        val missense = MutationsTest.mutation(false, Impact.MISSENSE)
        val victim = MutationsGene.tsgSampleSummary(DndsMutationTest.GENE, listOf(missense))
        MutationsTest.assertMutations(0, 1, victim.missense)
        MutationsTest.assertEmpty(victim.nonsense)
        MutationsTest.assertEmpty(victim.splice)
        MutationsTest.assertEmpty(victim.inframe)
        MutationsTest.assertEmpty(victim.frameshift)
    }

    @Test
    fun testMultiUnknownHit() {
        val missense = MutationsTest.mutation(false, Impact.MISSENSE)
        val victim = MutationsGene.tsgSampleSummary(DndsMutationTest.GENE, listOf(missense, missense))
        MutationsTest.assertRedundant(0, victim)
        MutationsTest.assertMutations(0, 2, victim.missense)
        MutationsTest.assertEmpty(victim.nonsense)
        MutationsTest.assertEmpty(victim.splice)
        MutationsTest.assertEmpty(victim.inframe)
        MutationsTest.assertEmpty(victim.frameshift)
    }

    @Test
    fun testSingleKnownHit() {
        val missense = MutationsTest.mutation(true, Impact.MISSENSE)
        val victim = MutationsGene.tsgSampleSummary(DndsMutationTest.GENE, listOf(missense))
        MutationsTest.assertRedundant(0, victim)
        MutationsTest.assertMutations(1, 0, victim.missense)
        MutationsTest.assertEmpty(victim.nonsense)
        MutationsTest.assertEmpty(victim.splice)
        MutationsTest.assertEmpty(victim.inframe)
        MutationsTest.assertEmpty(victim.frameshift)
    }

    @Test
    fun testMultiKnownHit() {
        val missense = MutationsTest.mutation(true, Impact.MISSENSE)
        val victim = MutationsGene.tsgSampleSummary(DndsMutationTest.GENE, listOf(missense, missense))
        MutationsTest.assertRedundant(1, victim)
        MutationsTest.assertMutations(1, 0, victim.missense)
        MutationsTest.assertEmpty(victim.nonsense)
        MutationsTest.assertEmpty(victim.splice)
        MutationsTest.assertEmpty(victim.inframe)
        MutationsTest.assertEmpty(victim.frameshift)
    }
}