package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.MM_PREFIX;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.MM_SUFFIX;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.getMMSkipValuesFromMMValue;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.getMMValueFromModCReadIndices;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.getModCReadIndicesFromMMSkipValues;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.getModCReadIndicesFromMMValue;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.SortedSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class BiomodalBamUtilsTest
{
    @Test
    public void testGetMMValueFromModCReadIndices()
    {
        // no modC
        String mmValue = getMMValueFromModCReadIndices("A".repeat(10).getBytes(), Sets.newTreeSet(), true);
        assertEquals(MM_PREFIX + MM_SUFFIX, mmValue);

        // single modC forward
        String readStr = "AAACCACA";
        SortedSet<Integer> modCReadIndices = Sets.newTreeSet();
        modCReadIndices.add(4);
        mmValue = getMMValueFromModCReadIndices(readStr.getBytes(), modCReadIndices, true);

        String expectedMMValue = MM_PREFIX + ",1" + MM_SUFFIX;
        assertEquals(expectedMMValue, mmValue);

        // single modC reverse
        readStr = "AAAGGAGA";
        modCReadIndices = Sets.newTreeSet();
        modCReadIndices.add(3);

        mmValue = getMMValueFromModCReadIndices(readStr.getBytes(), modCReadIndices, false);

        expectedMMValue = MM_PREFIX + ",2" + MM_SUFFIX;
        assertEquals(expectedMMValue, mmValue);

        // two modC forward
        readStr = "AAACCACA";
        modCReadIndices = Sets.newTreeSet();
        modCReadIndices.add(3);
        modCReadIndices.add(6);
        mmValue = getMMValueFromModCReadIndices(readStr.getBytes(), modCReadIndices, true);

        expectedMMValue = MM_PREFIX + ",0,1" + MM_SUFFIX;
        assertEquals(expectedMMValue, mmValue);

        // two modC reverse
        readStr = "AAAGGAGA";
        modCReadIndices = Sets.newTreeSet();
        modCReadIndices.add(4);
        modCReadIndices.add(6);
        mmValue = getMMValueFromModCReadIndices(readStr.getBytes(), modCReadIndices, false);

        expectedMMValue = MM_PREFIX + ",0,0" + MM_SUFFIX;
        assertEquals(expectedMMValue, mmValue);
    }

    @Test
    public void testGetModCReadIndicesFromMMSkipValues()
    {
        // no modC
        SortedSet<Integer> modCReadIndices = getModCReadIndicesFromMMSkipValues("A".repeat(10).getBytes(), Lists.newArrayList(), true);
        assertTrue(modCReadIndices.isEmpty());

        // single modC forward
        String readStr = "AAACCACA";
        List<Integer> skipValues = Lists.newArrayList(1);
        modCReadIndices = getModCReadIndicesFromMMSkipValues(readStr.getBytes(), skipValues, true);

        SortedSet<Integer> expectedModCReadIndices = Sets.newTreeSet();
        expectedModCReadIndices.add(4);
        assertEquals(expectedModCReadIndices, modCReadIndices);

        // single modC reverse
        readStr = "AAAGGAGA";
        skipValues = Lists.newArrayList(2);
        modCReadIndices = getModCReadIndicesFromMMSkipValues(readStr.getBytes(), skipValues, false);

        expectedModCReadIndices = Sets.newTreeSet();
        expectedModCReadIndices.add(3);
        assertEquals(expectedModCReadIndices, modCReadIndices);

        // two modC forward
        readStr = "AAACCACA";
        skipValues = Lists.newArrayList(0, 1);
        modCReadIndices = getModCReadIndicesFromMMSkipValues(readStr.getBytes(), skipValues, true);

        expectedModCReadIndices = Sets.newTreeSet();
        expectedModCReadIndices.add(3);
        expectedModCReadIndices.add(6);
        assertEquals(expectedModCReadIndices, modCReadIndices);

        // two modC reverse
        readStr = "AAAGGAGA";
        skipValues = Lists.newArrayList(0, 0);
        modCReadIndices = getModCReadIndicesFromMMSkipValues(readStr.getBytes(), skipValues, false);

        expectedModCReadIndices = Sets.newTreeSet();
        expectedModCReadIndices.add(4);
        expectedModCReadIndices.add(6);
        assertEquals(expectedModCReadIndices, modCReadIndices);
    }

    @Test
    public void testGetMMSkipValuesFromMMValue()
    {
        String mmValue = MM_PREFIX + MM_SUFFIX;
        List<Integer> skipValues = getMMSkipValuesFromMMValue(mmValue);
        assertTrue(skipValues.isEmpty());

        mmValue = MM_PREFIX + ",1" + MM_SUFFIX;
        skipValues = getMMSkipValuesFromMMValue(mmValue);
        List<Integer> expectedSkipValues = Lists.newArrayList(1);
        assertEquals(expectedSkipValues, skipValues);

        mmValue = MM_PREFIX + ",1,0" + MM_SUFFIX;
        skipValues = getMMSkipValuesFromMMValue(mmValue);
        expectedSkipValues = Lists.newArrayList(1, 0);
        assertEquals(expectedSkipValues, skipValues);
    }

    @Test
    public void testGetModCReadIndicesFromMMValue()
    {
        String readStr = "AAACCACA";
        String mmValue = MM_PREFIX + ",1,0" + MM_SUFFIX;
        SortedSet<Integer> modCReadIndices = getModCReadIndicesFromMMValue(readStr.getBytes(), mmValue, true);

        SortedSet<Integer> expectedModCReadIndices = Sets.newTreeSet();
        expectedModCReadIndices.add(4);
        expectedModCReadIndices.add(6);
        assertEquals(expectedModCReadIndices, modCReadIndices);
    }
}
