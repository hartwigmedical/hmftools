package com.hartwig.hmftools.paddle.mutation

import com.hartwig.hmftools.paddle.Impact
import com.hartwig.hmftools.paddle.dnds.DndsMutationTest.Companion.GENE
import com.hartwig.hmftools.paddle.mutation.MutationsTest.assertEmpty
import com.hartwig.hmftools.paddle.mutation.MutationsTest.assertMutations
import com.hartwig.hmftools.paddle.mutation.MutationsTest.assertRedundant
import com.hartwig.hmftools.paddle.mutation.MutationsTest.mutation
import com.hartwig.hmftools.paddle.mutation.MutationsTest.synonymousMutation
import org.junit.Test

class OncoMutationTest {

    @Test
    fun testSingleUnknownMissense() {
        val missense = mutation(false, Impact.MISSENSE)
        val victim = MutationsGene.oncoSampleSummary(GENE, listOf(missense))
        assertMutations(0, 1, victim.missense)
        assertEmpty(victim.nonsense)
        assertEmpty(victim.splice)
        assertEmpty(victim.inframe)
        assertEmpty(victim.frameshift)
    }

    @Test
    fun testMultiUnknownMissense() {
        val missense = mutation(false, Impact.MISSENSE)
        val victim = MutationsGene.oncoSampleSummary(GENE, listOf(missense, missense))
        assertRedundant(1, victim)
        assertMutations(0, 1, victim.missense)
        assertEmpty(victim.nonsense)
        assertEmpty(victim.splice)
        assertEmpty(victim.inframe)
        assertEmpty(victim.frameshift)
    }

    @Test
    fun testSingleKnownMissense() {
        val missense = mutation(true, Impact.MISSENSE)
        val victim = MutationsGene.oncoSampleSummary(GENE, listOf(missense))
        assertMutations(1, 0, victim.missense)
        assertEmpty(victim.nonsense)
        assertEmpty(victim.splice)
        assertEmpty(victim.inframe)
        assertEmpty(victim.frameshift)
    }

    @Test
    fun testMultipleMissense() {
        val missense = mutation(true, Impact.MISSENSE)
        val victim = MutationsGene.oncoSampleSummary(GENE, listOf(missense, missense))
        assertRedundant(1, victim)
        assertMutations(1, 0, victim.missense)
        assertEmpty(victim.nonsense)
        assertEmpty(victim.splice)
        assertEmpty(victim.inframe)
        assertEmpty(victim.frameshift)
    }

    @Test
    fun testDifferentTypes() {
        val missense = mutation(false, Impact.MISSENSE)
        val nonsense = mutation(true, Impact.NONSENSE)
        val synonymous = synonymousMutation()
        val victim = MutationsGene.oncoSampleSummary(GENE, listOf(missense, nonsense, synonymous))
        assertRedundant(1, victim)
        assertMutations(1, 0, victim.nonsense)
        assertEmpty(victim.missense)
        assertEmpty(victim.splice)
        assertEmpty(victim.inframe)
        assertEmpty(victim.frameshift)
    }
}