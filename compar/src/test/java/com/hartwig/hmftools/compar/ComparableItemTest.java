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

    // Builder for objects to be tested. Alternate initializer in builder needs to contain a value different from the default for every
    // field that is compared in the "findMismatch" method, used in the "matches" method or that is used in the "reportable" method.
    protected TestComparableItemBuilder<B, I> builder;

    // Map from name of compared field in "findMismatch" to initializer that causes a mismatch in that field.
    protected Map<String, Consumer<B>> fieldToAlternateValueInitializer;

    // Map from name/description of field used in "matches" method to initializer that causes a mismatch in that field
    protected Map<String, Consumer<B>> nameToAlternateIndexInitializer;

    // Map from name/description of field used in "reportable" method to initializer that causes a mismatch in that field
    protected Map<String, Consumer<B>> reportabilityFieldToFalseReportabilityInitializer;

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
