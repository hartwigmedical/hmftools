package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.ComparTestUtil.assertDifferencesAreForFields;
import static com.hartwig.hmftools.compar.ComparTestUtil.assertSingleFieldMismatch;
import static com.hartwig.hmftools.compar.ComparTestUtil.assertValueDifferencesAsExpected;
import static com.hartwig.hmftools.compar.ComparTestUtil.combine;
import static com.hartwig.hmftools.compar.ComparTestUtil.union;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Test;

public abstract class ComparableItemTest<I extends ComparableItem, C extends ItemComparer, B>
{
    protected C comparer;

    // Builder for objects to be tested. Default object needs to be "reportable" and "pass".
    // Alternate initializer in builder needs to contain a value different from the default for every
    // field that is compared in the "findMismatch" method, used in the "matches" method or that is used in the "reportable" method.
    // Alternate object should still be pass.
    protected TestComparableItemBuilder<B, I> builder;

    // Map from name of field compared in "findMismatch" to initializer that changes that field from the default value.
    // This is meant for causing "VALUE" differences in testing.
    protected Map<String, Consumer<B>> fieldToAlternateValueInitializer;

    // Map from name of field or fields used in "matches" method to initializer that changes that field from the default value.
    // This is meant for creating non-matching objects in testing.
    protected Map<String, Consumer<B>> nameToAlternateIndexInitializer;

    // Map from name of field used in "reportable" method to initializer that changes that field from the default value.
    // This is meant for causing differences in reportability in testing.
    protected Map<String, Consumer<B>> reportabilityFieldToFalseReportabilityInitializer;

    // Map from name of field or fields used in "isPass" method to initializer that causes "isPass" to be false.
    // This is meant for causing differences in PASS status in testing.
    protected Map<String, Consumer<B>> nameToNonPassInitializer;

    @Test
    public void fullyMatchesSelfInDetailedMode()
    {
        assertFullyMatchesSelf(MatchLevel.DETAILED);
    }

    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        assertFullyMatchesSelf(MatchLevel.REPORTABLE);
    }

    @Test
    public void fullyDifferent()
    {
        assertValueDifferencesAsExpected(
                builder.create(),
                builder.createWithAlternateDefaults(),
                MatchLevel.DETAILED,
                createDefaultThresholds(),
                getAllValueFieldNames(),
                nameToAlternateIndexInitializer.isEmpty()
        );
    }

    @Test
    public void singleIndexMismatchesAreRecognized()
    {
        I victim = builder.create();
        for(Map.Entry<String, Consumer<B>> entry : nameToAlternateIndexInitializer.entrySet())
        {
            I indexMismatch = builder.create(entry.getValue());
            assertFalse("Test mismatch in index field (victim.matches(alt)): " + entry.getKey(), victim.matches(indexMismatch));
            assertFalse("Test mismatch in index field (alt.matches(victim)): " + entry.getKey(), indexMismatch.matches(victim));
        }
    }

    @Test
    public void singleFieldMismatchesAreRecognizedInDetailedMode()
    {
        assertSingleFieldMismatchesAreRecognized(MatchLevel.DETAILED);
    }

    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        assertSingleFieldMismatchesAreRecognized(MatchLevel.REPORTABLE);
    }

    @Test
    public void reportabilityMismatchesAreRecognized()
    {
        for(Map.Entry<String, Consumer<B>> entry : reportabilityFieldToFalseReportabilityInitializer.entrySet())
        {
            String field = entry.getKey();
            I reportableVictim = builder.create();
            I nonReportableVictim = builder.create(entry.getValue());
            DiffThresholds diffThresholds = createDefaultThresholds();

            assertSingleFieldMismatch(field, reportableVictim, nonReportableVictim, MatchLevel.DETAILED, diffThresholds, MismatchType.VALUE);
            assertSingleFieldMismatch(field, nonReportableVictim, reportableVictim, MatchLevel.DETAILED, diffThresholds, MismatchType.VALUE);

            assertSingleFieldMismatch(field, reportableVictim, nonReportableVictim, MatchLevel.REPORTABLE, diffThresholds, MismatchType.REF_ONLY);
            assertSingleFieldMismatch(field, nonReportableVictim, reportableVictim, MatchLevel.REPORTABLE, diffThresholds, MismatchType.NEW_ONLY);
        }
    }

    @Test
    public void onlyMatchesIndex()
    {
        DiffThresholds diffThresholds = createDefaultThresholds();

        I refVictim = builder.create();
        I newVictim = builder.create(combine(getAllAlternateValueInitializers()));

        assertTrue(refVictim.matches(newVictim));
        Mismatch mismatch = refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.Type);
        assertEquals(refVictim, mismatch.RefItem);
        assertEquals(newVictim, mismatch.NewItem);
        assertDifferencesAreForFields(getAllValueFieldNames(), mismatch.DiffValues);
    }

    @Test
    public void categoriesMatch()
    {
        assertEquals(comparer.category(), builder.create().category());
    }

    @Test
    public void amountOfDisplayValuesMatches()
    {
        assertEquals(comparer.comparedFieldNames().size(), builder.create().displayValues().size());
    }

    @Test
    public void hasKeyIfItShould()
    {
        assertEquals(nameToAlternateIndexInitializer.isEmpty(), builder.create().key().isEmpty());
    }

    @Test
    public void doubleNonPassIsIgnored()
    {
        for(Map.Entry<String, Consumer<B>> entry : nameToNonPassInitializer.entrySet())
        {
            String name = entry.getKey();
            final Consumer<B> initializer = entry.getValue();
            I refVictim = builder.create(initializer);
            I newVictim = builder.createWithAlternateDefaults(initializer);

            DiffThresholds diffThresholds = createDefaultThresholds();
            assertNull("Test non-PASS due to " + name + " is ignored when not including matches",
                    refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, false));
            assertNull("Test non-PASS due to " + name + " is ignored when including matches",
                    refVictim.findMismatch(newVictim, MatchLevel.DETAILED, diffThresholds, true));
        }
    }

    @Test
    public void nonPassIsNotCallInDetailedMode()
    {
        for(Map.Entry<String, Consumer<B>> entry : nameToNonPassInitializer.entrySet())
        {
            String name = entry.getKey();
            final Consumer<B> initializer = entry.getValue();
            I passVictim = builder.create();
            I nonPassVictim = builder.createWithAlternateDefaults(initializer);

            DiffThresholds diffThresholds = createDefaultThresholds();
            assertEquals("Test non-PASS due to " + name + " can cause REF_ONLY when not including matches",
                    MismatchType.REF_ONLY, passVictim.findMismatch(nonPassVictim, MatchLevel.DETAILED, diffThresholds, false).Type);
            assertEquals("Test non-PASS due to " + name + " can cause REF_ONLY when including matches",
                    MismatchType.REF_ONLY, passVictim.findMismatch(nonPassVictim, MatchLevel.DETAILED, diffThresholds, true).Type);

            assertEquals("Test non-PASS due to " + name + " can cause NEW_ONLY when not including matches",
                    MismatchType.NEW_ONLY, nonPassVictim.findMismatch(passVictim, MatchLevel.DETAILED, diffThresholds, false).Type);
            assertEquals("Test non-PASS due to " + name + " can cause NEW_ONLY when including matches",
                    MismatchType.NEW_ONLY, nonPassVictim.findMismatch(passVictim, MatchLevel.DETAILED, diffThresholds, true).Type);
        }
    }

    private void assertFullyMatchesSelf(final MatchLevel matchLevel)
    {
        DiffThresholds diffThresholds = createDefaultThresholds();

        I victim = builder.create();

        assertTrue(victim.matches(victim));

        assertNull(victim.findMismatch(victim, matchLevel, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, victim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(victim, matchLevel, diffThresholds, true));
    }

    private void assertSingleFieldMismatchesAreRecognized(final MatchLevel matchLevel)
    {
        for(Map.Entry<String, Consumer<B>> entry : fieldToAlternateValueInitializer.entrySet())
        {
            String field = entry.getKey();
            I refVictim = builder.create();
            I newVictim = builder.create(entry.getValue());
            DiffThresholds diffThresholds = createDefaultThresholds();

            assertSingleFieldMismatch(field, refVictim, newVictim, matchLevel, diffThresholds, MismatchType.VALUE);
        }
    }

    protected DiffThresholds createDefaultThresholds()
    {
        DiffThresholds diffThresholds = new DiffThresholds();
        comparer.registerThresholds(diffThresholds);
        return diffThresholds;
    }

    private Set<String> getAllValueFieldNames()
    {
        return union(fieldToAlternateValueInitializer.keySet(), reportabilityFieldToFalseReportabilityInitializer.keySet());
    }

    private Set<Consumer<B>> getAllAlternateValueInitializers()
    {
        return union(fieldToAlternateValueInitializer.values(), reportabilityFieldToFalseReportabilityInitializer.values());
    }
}
