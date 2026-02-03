package com.hartwig.hmftools.cider

import kotlin.test.*

class SomaticHypermutationTest {
    @Test
    fun testDetermineFromVGeneIdentity()
    {
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(98.0))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(98.01))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(99.0))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(99.99))
        assertEquals(SomaticHypermutationStatus.UNMUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(100.0))

        assertEquals(SomaticHypermutationStatus.MUTATED_BORDERLINE, SomaticHypermutationStatus.determineFromVGeneIdentity(97.0))
        assertEquals(SomaticHypermutationStatus.MUTATED_BORDERLINE, SomaticHypermutationStatus.determineFromVGeneIdentity(97.01))
        assertEquals(SomaticHypermutationStatus.MUTATED_BORDERLINE, SomaticHypermutationStatus.determineFromVGeneIdentity(97.99))

        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(96.99))
        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(96.0))
        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(95.0))
        assertEquals(SomaticHypermutationStatus.MUTATED, SomaticHypermutationStatus.determineFromVGeneIdentity(90.0))
    }
}