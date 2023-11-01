package com.hartwig.hmftools.dnds.calcs;

public class TsgMutationTest
{
    /*
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

     */

    /*

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
     */
}
