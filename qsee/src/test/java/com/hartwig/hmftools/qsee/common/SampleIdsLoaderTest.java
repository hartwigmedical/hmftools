package com.hartwig.hmftools.qsee.common;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class SampleIdsLoaderTest
{
    @Test
    public void canLoadFromStringsWithMatchingCounts()
    {
        SampleIdsLoader loader = new SampleIdsLoader();
        loader.fromStrings("TUMOR001,TUMOR002", "REF001,REF002");

        assertEquals(2, loader.tumorIds().size());
        assertEquals(2, loader.referenceIds().size());
    }

    @Test
    public void canLoadFromStringsWithTumorOnly()
    {
        SampleIdsLoader loader = new SampleIdsLoader();
        loader.fromStrings("TUMOR001,TUMOR002", null);

        assertEquals(2, loader.tumorIds().size());
        assertEquals(0, loader.referenceIds().size());
    }

    @Test(expected = IllegalStateException.class)
    public void cannotLoadFromStringsWithMismatchedCounts()
    {
        SampleIdsLoader loader = new SampleIdsLoader();
        loader.fromStrings("TUMOR001,TUMOR002", "REF001");
    }

    @Test
    public void canLoadFromLinesWithMatchingCounts()
    {
        List<String> lines = new ArrayList<>();
        lines.add("TumorId\tReferenceId");
        lines.add("TUMOR001\tREF001");
        lines.add("TUMOR002\tREF002");

        SampleIdsLoader loader = new SampleIdsLoader();
        loader.fromLines(lines);

        assertEquals(2, loader.tumorIds().size());
        assertEquals(2, loader.referenceIds().size());
    }

    @Test
    public void canLoadFromLinesWithTumorOnly()
    {
        List<String> lines = new ArrayList<>();
        lines.add("TumorId");
        lines.add("TUMOR001");
        lines.add("TUMOR002");

        SampleIdsLoader loader = new SampleIdsLoader();
        loader.fromLines(lines);

        assertEquals(2, loader.tumorIds().size());
        assertEquals(0, loader.referenceIds().size());
    }

    @Test(expected = IllegalStateException.class)
    public void cannotLoadFromLinesWithMismatchedCounts()
    {
        List<String> lines = new ArrayList<>();
        lines.add("TumorId\tReferenceId");
        lines.add("TUMOR001\tREF001");
        lines.add("TUMOR002\t");

        SampleIdsLoader loader = new SampleIdsLoader();
        loader.fromLines(lines);
    }
}