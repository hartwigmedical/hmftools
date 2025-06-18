package com.hartwig.hmftools.compar.mutation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Test;

public class SomaticVariantDataTest
{
    @Test
    public void fullyMatchesSelf()
    {
        var victim = SomaticVariantDataTestFactory.createDefault().build();
        var diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(victim));
        assertNull(victim.findMismatch(victim, MatchLevel.DETAILED, diffThresholds, false));
        assertNull(victim.findMismatch(victim, MatchLevel.REPORTABLE, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, victim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(victim, MatchLevel.DETAILED, diffThresholds, true));
        assertEquals(expectedMatch, victim.findMismatch(victim, MatchLevel.REPORTABLE, diffThresholds, true));
    }

    @Test
    public void fullyMatchesSelfWithLiftover()
    {
        var victim = SomaticVariantDataTestFactory.createDefault().withComparisonChromosome("8").withComparisonPosition(10000).build();
        var liftoverVictim = SomaticVariantDataTestFactory.createDefault().withChromosome("8").withPosition(10000).build();
        var diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.REPORTABLE, diffThresholds, true));
    }

    @Test
    public void onlyMatchesIndex()
    {
        var refVictim = SomaticVariantDataTestFactory.createDefault().build();
        var newVictim = SomaticVariantDataTestFactory.createAlternateDefault().withChromosome(refVictim.Chromosome)
                .withPosition(refVictim.Position)
                .withRef(refVictim.Ref)
                .withAlt(refVictim.Alt)
                .withType(refVictim.Type)
                .withComparisonChromosome(refVictim.mComparisonChromosome)
                .withComparisonPosition(refVictim.mComparisonPosition)
                .build();

        var diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));
        var mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        assertEquals(18, mismatch.DiffValues().size());
    }

    @Test
    public void nonPurpleMatchHandledCorrectly()
    {
        var refVictim = SomaticVariantDataTestFactory.createDefault().withHasPurpleAnnotation(false).build();
        var newVictim = SomaticVariantDataTestFactory.createAlternateDefault().withChromosome(refVictim.Chromosome)
                .withPosition(refVictim.Position)
                .withRef(refVictim.Ref)
                .withAlt(refVictim.Alt)
                .withType(refVictim.Type)
                .withComparisonChromosome(refVictim.mComparisonChromosome)
                .withComparisonPosition(refVictim.mComparisonPosition)
                .withHasPurpleAnnotation(false)
                .build();

        var diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));
        var mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        assertEquals(12, mismatch.DiffValues().size());
    }

    @Test
    public void unfilteredMatchHandledCorrectly()
    {
        var passVictim = SomaticVariantDataTestFactory.createDefault().build();
        var filteredVictim = SomaticVariantDataTestFactory.createAlternateDefault().withChromosome(passVictim.Chromosome)
                .withPosition(passVictim.Position)
                .withRef(passVictim.Ref)
                .withAlt(passVictim.Alt)
                .withType(passVictim.Type)
                .withComparisonChromosome(passVictim.mComparisonChromosome)
                .withComparisonPosition(passVictim.mComparisonPosition)
                .withIsFromUnfilteredVcf(true)
                .withHasPurpleAnnotation(false)
                .build();

        var diffThresholds = createDefaultThresholds();

        assertTrue(passVictim.matches(filteredVictim));
        var mismatch = passVictim.findMismatch(filteredVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.REF_ONLY, mismatch.MismatchType());
        assertEquals(passVictim, mismatch.RefItem());
        assertEquals(filteredVictim, mismatch.NewItem());
        assertEquals(7, mismatch.DiffValues().size());

        assertTrue(filteredVictim.matches(passVictim));
        var oppositeMismatch = filteredVictim.findMismatch(passVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.NEW_ONLY, oppositeMismatch.MismatchType());
        assertEquals(filteredVictim, oppositeMismatch.RefItem());
        assertEquals(passVictim, oppositeMismatch.NewItem());
        assertEquals(7, oppositeMismatch.DiffValues().size());
    }

    @Test
    public void indexMismatchesAreRecognized()
    {
        var baseVictim = SomaticVariantDataTestFactory.createDefault();
        var victim = baseVictim.build();

        var alternateVictim = SomaticVariantDataTestFactory.createAlternateDefault().build();

        var chromosomeMismatch =
                baseVictim.withChromosome(alternateVictim.Chromosome).withComparisonChromosome(alternateVictim.Chromosome).build();

        assertFalse(victim.matches(chromosomeMismatch));
        assertFalse(chromosomeMismatch.matches(victim));

        var positionMismatch = baseVictim.withPosition(alternateVictim.Position).withComparisonPosition(alternateVictim.Position).build();

        assertFalse(victim.matches(positionMismatch));
        assertFalse(positionMismatch.matches(victim));

        var refMismatch = baseVictim.withRef(alternateVictim.Ref).build();

        assertFalse(victim.matches(refMismatch));
        assertFalse(refMismatch.matches(victim));

        var altMismatch = baseVictim.withAlt(alternateVictim.Alt).build();

        assertFalse(victim.matches(altMismatch));
        assertFalse(altMismatch.matches(victim));

        var variantTypeMismatch = baseVictim.withType(alternateVictim.Type).build();

        assertFalse(victim.matches(variantTypeMismatch));
        assertFalse(variantTypeMismatch.matches(victim));
    }

    private static DiffThresholds createDefaultThresholds()
    {
        var config = new ComparConfig();
        var comparer = new SomaticVariantComparer(config);

        var diffThresholds = new DiffThresholds();
        comparer.registerThresholds(diffThresholds);
        return diffThresholds;
    }
}
