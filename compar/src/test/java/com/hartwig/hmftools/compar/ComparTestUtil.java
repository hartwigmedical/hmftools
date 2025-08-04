package com.hartwig.hmftools.compar;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

public class ComparTestUtil
{
    public static void assertDifferencesAreForFields(final Set<String> expectedFields, final List<String> differences)
    {
        List<String> sortedExpectedFields = expectedFields.stream().sorted().toList();
        List<String> sortedObservedFields = differences.stream()
                .map(ComparTestUtil::extractFieldNameFromDifference)
                .distinct()
                .sorted()
                .toList();
        assertEquals(sortedExpectedFields, sortedObservedFields);
    }

    public static void assertSingleFieldMismatch(final String field, final ComparableItem refVictim, final ComparableItem newVictim,
            final MatchLevel matchLevel, final DiffThresholds diffThresholds, final MismatchType expectedMismatchType)
    {
        String testMessage = "Test mismatch in field '" + field + "' at match level '" + matchLevel + "'";
        assertTrue(testMessage, refVictim.matches(newVictim));
        assertTrue(testMessage, newVictim.matches(refVictim));

        Mismatch detailedMismatch = refVictim.findMismatch(newVictim, matchLevel, diffThresholds, false);

        assertEquals(testMessage, expectedMismatchType, detailedMismatch.MismatchType());
        assertEquals(testMessage, refVictim, detailedMismatch.RefItem());
        assertEquals(testMessage, newVictim, detailedMismatch.NewItem());
        assertEquals(testMessage, 1, detailedMismatch.DiffValues().size());
        assertEquals(testMessage, field, extractFieldNameFromDifference(detailedMismatch.DiffValues().get(0)));
    }

    public static void assertValueDifferencesAsExpected(final ComparableItem refVictim, final ComparableItem newVictim,
            final MatchLevel matchLevel, final DiffThresholds diffThresholds, final Set<String> expectedFieldNames,
            final boolean expectIndexMatch)
    {
        if(expectIndexMatch)
        {
            assertTrue(refVictim.matches(newVictim));
        }
        else
        {
            assertFalse(refVictim.matches(newVictim));
        }

        Mismatch mismatch = refVictim.findMismatch(newVictim, matchLevel, diffThresholds, false);

        assertEquals(MismatchType.VALUE, mismatch.MismatchType());
        assertEquals(refVictim, mismatch.RefItem());
        assertEquals(newVictim, mismatch.NewItem());

        assertDifferencesAreForFields(expectedFieldNames, mismatch.DiffValues());
    }

    public static String extractFieldNameFromDifference(final String difference)
    {
        return difference.split("\\(")[0];
    }

    public static <T> Set<T> union(final List<Collection<T>> sets)
    {
        return sets.stream().flatMap(Collection::stream).collect(Collectors.toSet());
    }

    public static <T> Set<T> union(final Collection<T> set1, final Collection<T> set2)
    {
        return union(List.of(set1, set2));
    }

    public static <T> Set<T> union(final Collection<T> set1, final Collection<T> set2, final Collection<T> set3)
    {
        return union(List.of(set1, set2, set3));
    }

    public static <T> Consumer<T> combine(final Collection<Consumer<T>> initializers)
    {
        return initializers.stream().reduce(s -> {}, Consumer::andThen);
    }

    public static <T> Consumer<T> combine(final Consumer<T> initializer1, final Consumer<T> initializer2)
    {
        return combine(List.of(initializer1, initializer2));
    }
}
