package com.hartwig.hmftools.compar.linx;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;
import java.util.function.Consumer;

import org.junit.Test;

public class BreakendDataTest
{
    @Test
    public void matchesSelf()
    {
        BreakendData victim = TestBreakendDataBuilder.BUILDER.create();
        assertTrue(victim.matches(victim));
    }

    @Test
    public void doesNotMatchOther()
    {
        BreakendData victim = TestBreakendDataBuilder.BUILDER.create();
        BreakendData otherVictim = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();
        assertFalse(victim.matches(otherVictim));
        assertFalse(otherVictim.matches(victim));
    }

    @Test
    public void canMatchAfterLiftover()
    {
        BreakendData refVictim = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.comparisonChromosome = "21";
            b.comparisonPosition = 40000000;
        });
        BreakendData newVictim = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.chromosome = "21";
            b.position = 40000000;
        });
        assertTrue(refVictim.matches(newVictim));
        assertTrue(newVictim.matches(refVictim));
    }

    @Test
    public void transcriptIsCheckedWhenAtLeastOneNonCanonical()
    {
        BreakendData victim = TestBreakendDataBuilder.BUILDER.create();
        BreakendData victimWithNonCanonical = TestBreakendDataBuilder.BUILDER.create(b -> b.canonical = false);
        BreakendData victimWithOtherCanonical = TestBreakendDataBuilder.BUILDER.create(b -> b.transcriptId = "OTHER");
        BreakendData victimWithOtherNonCanonical = TestBreakendDataBuilder.BUILDER.create(b -> {
            b.transcriptId = "OTHER";
            b.canonical = false;
        });

        assertTrue(victim.matches(victimWithNonCanonical));
        assertTrue(victimWithNonCanonical.matches(victim));

        assertFalse(victim.matches(victimWithOtherNonCanonical));
        assertFalse(victimWithOtherNonCanonical.matches(victim));

        assertTrue(victimWithOtherCanonical.matches(victimWithOtherNonCanonical));
        assertTrue(victimWithOtherNonCanonical.matches(victimWithOtherCanonical));

        assertFalse(victimWithNonCanonical.matches(victimWithOtherNonCanonical));
        assertFalse(victimWithOtherNonCanonical.matches(victimWithNonCanonical));
    }

    @Test
    public void transcriptIsNotCheckedWhenOnlyCanonical()
    {
        BreakendData victim = TestBreakendDataBuilder.BUILDER.create();
        BreakendData victimWithOtherCanonical = TestBreakendDataBuilder.BUILDER.create(b -> b.transcriptId = "OTHER");

        assertTrue(victim.matches(victimWithOtherCanonical));
        assertTrue(victimWithOtherCanonical.matches(victim));
    }

    @Test
    public void breakendHomologyIsRespected()
    {
        BreakendData victim = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();

        int[] otherHomologyOffset = new int[] { -5, 4 };
        BreakendData victimBarelyOverlapFromAbove = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults(b -> {
            b.position = victim.ComparisonPosition + victim.HomologyOffset[1] - otherHomologyOffset[0];
            b.comparisonPosition = victim.Position + victim.HomologyOffset[1] - otherHomologyOffset[0];
            b.homologyOffset = otherHomologyOffset;
        });
        BreakendData victimBarelyNoOverlapFromAbove = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults(b -> {
            b.position = victim.ComparisonPosition + victim.HomologyOffset[1] - otherHomologyOffset[0] + 1;
            b.comparisonPosition = victim.Position + victim.HomologyOffset[1] - otherHomologyOffset[0] + 1;
            b.homologyOffset = otherHomologyOffset;
        });
        BreakendData victimBarelyOverlapFromBelow = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults(b -> {
            b.position = victim.ComparisonPosition + victim.HomologyOffset[0] - otherHomologyOffset[1];
            b.comparisonPosition = victim.Position + victim.HomologyOffset[0] - otherHomologyOffset[1];
            b.homologyOffset = otherHomologyOffset;
        });
        BreakendData victimBarelyNoOverlapFromBelow = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults(b -> {
            b.position = victim.ComparisonPosition + victim.HomologyOffset[0] - otherHomologyOffset[1] - 1;
            b.comparisonPosition = victim.Position + victim.HomologyOffset[0] - otherHomologyOffset[1] - 1;
            b.homologyOffset = otherHomologyOffset;
        });

        assertTrue(victim.matches(victim));

        assertTrue(victim.matches(victimBarelyOverlapFromAbove));
        assertTrue(victimBarelyOverlapFromAbove.matches(victim));

        assertTrue(victim.matches(victimBarelyOverlapFromBelow));
        assertTrue(victimBarelyOverlapFromBelow.matches(victim));

        assertFalse(victim.matches(victimBarelyNoOverlapFromAbove));
        assertFalse(victimBarelyNoOverlapFromAbove.matches(victim));

        assertFalse(victim.matches(victimBarelyNoOverlapFromBelow));
        assertFalse(victimBarelyNoOverlapFromBelow.matches(victim));
    }

    @Test
    public void singleIndexMismatchesInBreakendAreRecognized()
    {
        BreakendData alternateDefaults = TestBreakendDataBuilder.BUILDER.createWithAlternateDefaults();

        Map<String, Consumer<TestBreakendDataBuilder>> fieldToAlternateValueInitializer = Map.of(
                "SvType", b -> b.svType = alternateDefaults.SvType,
                "Chromosome", b ->
                {
                    b.chromosome = alternateDefaults.Chromosome;
                    b.comparisonChromosome = alternateDefaults.ComparisonChromosome;
                },
                "Orientation", b -> b.orientation = alternateDefaults.Orientation,
                "Position", b -> {
                    b.position = alternateDefaults.Position;
                    b.comparisonPosition = alternateDefaults.ComparisonPosition;
                }
        );

        for(Map.Entry<String, Consumer<TestBreakendDataBuilder>> entry : fieldToAlternateValueInitializer.entrySet())
        {
            String field = entry.getKey();
            BreakendData defaultVictim = TestBreakendDataBuilder.BUILDER.create();
            BreakendData nonDefaultVictim = TestBreakendDataBuilder.BUILDER.create(entry.getValue());

            assertFalse("Test difference in " + field, defaultVictim.matches(nonDefaultVictim));
            assertFalse("Test difference in " + field, nonDefaultVictim.matches(defaultVictim));
        }
    }
}
