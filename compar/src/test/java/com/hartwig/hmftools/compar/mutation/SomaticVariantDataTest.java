package com.hartwig.hmftools.compar.mutation;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Set;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
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
        var victim = createVictim();
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
        var victim = createVictim().withComparisonChromosome("8").withComparisonPosition(10000);
        var liftoverVictim = createVictim().withChromosome("8").withPosition(10000);
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
        var refVictim = createVictim();
        var newVictim = createAlternateVictim().withChromosome(refVictim.Chromosome)
                .withPosition(refVictim.Position)
                .withRef(refVictim.Ref)
                .withAlt(refVictim.Alt)
                .withType(refVictim.Type)
                .withComparisonChromosome(refVictim.mComparisonChromosome)
                .withComparisonPosition(refVictim.mComparisonPosition);

        var diffThresholds = createDefaultThresholds();

        assertTrue(refVictim.matches(newVictim));
        var mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());
        assertEquals(18, mismatch.DiffValues().size());
    }

    @Test
    public void unfilteredMatchHandledCorrectly()
    {
        var passVictim = createVictim();
        var filteredVictim = createAlternateVictim().withChromosome(passVictim.Chromosome)
                .withPosition(passVictim.Position)
                .withRef(passVictim.Ref)
                .withAlt(passVictim.Alt)
                .withType(passVictim.Type)
                .withComparisonChromosome(passVictim.mComparisonChromosome)
                .withComparisonPosition(passVictim.mComparisonPosition)
                .withIsFromUnfilteredVcf(true);

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
        var victim = createVictim();
        var alternateVictim = createAlternateVictim();

        var chromosomeMismatch = victim.withChromosome(alternateVictim.Chromosome).withComparisonChromosome(alternateVictim.Chromosome);

        assertFalse(victim.matches(chromosomeMismatch));
        assertFalse(chromosomeMismatch.matches(victim));

        var positionMismatch = victim.withPosition(alternateVictim.Position).withComparisonPosition(alternateVictim.Position);

        assertFalse(victim.matches(positionMismatch));
        assertFalse(positionMismatch.matches(victim));

        var refMismatch = victim.withRef(alternateVictim.Ref);

        assertFalse(victim.matches(refMismatch));
        assertFalse(refMismatch.matches(victim));

        var altMismatch = victim.withAlt(alternateVictim.Alt);

        assertFalse(victim.matches(altMismatch));
        assertFalse(altMismatch.matches(victim));

        var variantTypeMismatch = victim.withType(alternateVictim.Type);

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

    private static SomaticVariantData createVictim()
    {
        return new SomaticVariantData("7", 140453136, "A", "T", VariantType.SNP, "BRAF", true,
                Hotspot.HOTSPOT, VariantTier.HOTSPOT, false, "missense_variant", "MISSENSE",
                "c.1799T>A", "p.Val600Glu", null, false,
                275, 0., Set.of("PASS"), 1.1, 0.45,
                new AllelicDepth(116, 21), false, "7", 140453136);
    }

    private static SomaticVariantData createAlternateVictim()
    {
        return new SomaticVariantData("8", 10000, "C", "G", VariantType.INDEL, "BRCA1", false,
                Hotspot.NEAR_HOTSPOT, VariantTier.PANEL, true, "synonymous_variant", "SYNONYMOUS",
                "c.1800T>A", "p.Val601Glu", "OTHER_EFFECT", true,
                512, 1., Set.of("PON"), 3.6, 1.1,
                new AllelicDepth(312, 50), false, "8", 10000);
    }
}
