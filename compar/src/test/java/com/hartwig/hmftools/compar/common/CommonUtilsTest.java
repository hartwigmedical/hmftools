package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.compar.ComparConfig.NEW_SOURCE;
import static com.hartwig.hmftools.compar.ComparConfig.REF_SOURCE;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonGenomePosition;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.ComparableItem;

import org.junit.Test;

public class CommonUtilsTest
{
    @Test
    public void emptyComparison()
    {
        DiffThresholds diffThresholds = new DiffThresholds();

        List<Mismatch> mismatches = new ArrayList<>();
        List<ComparableItem> refItems = new ArrayList<>();
        List<ComparableItem> newItems = new ArrayList<>();
        MatchLevel matchLevel = MatchLevel.DETAILED;
        boolean includeMatches = false;

        CommonUtils.compareItems(mismatches, matchLevel, diffThresholds, includeMatches, refItems, newItems);

        assertTrue(mismatches.isEmpty());
    }

    @Test
    public void testDetermineComparisonGenomePositionWithoutLiftover()
    {
        GenomeLiftoverCache lifoverCache = new GenomeLiftoverCache(true);

        assertEquals(new BasePosition("8", 10000),
                determineComparisonGenomePosition("8", 10000, REF_SOURCE, false, lifoverCache));
        assertEquals(new BasePosition("chr10", 10000),
                determineComparisonGenomePosition("chr10", 10000, REF_SOURCE, false, lifoverCache));
        assertEquals(new BasePosition("X", 10000),
                determineComparisonGenomePosition("X", 10000, NEW_SOURCE, false, lifoverCache));
        assertEquals(new BasePosition("chr3", 10000),
                determineComparisonGenomePosition("chr3", 10000, NEW_SOURCE, false, lifoverCache));
    }

    @Test
    public void testDetermineComparisonGenomePositionWithLiftover()
    {
        GenomeLiftoverCache lifoverCache = new GenomeLiftoverCache(true);

        assertTrue(lifoverCache.hasMappings());
        assertEquals(new BasePosition("chr8", 150000),
                determineComparisonGenomePosition("8", 100000, REF_SOURCE, true, lifoverCache));
        assertEquals(new BasePosition("chr10", 100000),
                determineComparisonGenomePosition("chr10", 100000, REF_SOURCE, true, lifoverCache));
        assertEquals(new BasePosition("X", 100000),
                determineComparisonGenomePosition("X", 100000, NEW_SOURCE, true, lifoverCache));
        assertEquals(new BasePosition("3", 141683),
                determineComparisonGenomePosition("chr3", 100000, NEW_SOURCE, true, lifoverCache));
    }

    @Test
    public void testDetermineComparisonChromosomeWithoutLiftover()
    {
        assertEquals("8", determineComparisonChromosome("8", false));
        assertEquals("chr10", determineComparisonChromosome("chr10", false));
        assertEquals("X", determineComparisonChromosome("X", false));
        assertEquals("chrY", determineComparisonChromosome("chrY", false));
  }

    @Test
    public void testDetermineComparisonChromosomeWithLiftover()
    {
        assertEquals("8", determineComparisonChromosome("8", true));
        assertEquals("10", determineComparisonChromosome("chr10", true));
        assertEquals("X", determineComparisonChromosome("X", true));
        assertEquals("Y", determineComparisonChromosome("chrY", true));
    }
}
