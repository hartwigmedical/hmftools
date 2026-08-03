package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.FieldCheckCache;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class NovelSpliceJunctionDataTest
        extends ComparableItemTest<NovelSpliceJunctionData, NovelSpliceJunctionComparer, TestNovelSpliceJunctionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new NovelSpliceJunctionComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestNovelSpliceJunctionDataBuilder.BUILDER;
        NovelSpliceJunctionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                NovelSpliceJunctionComparer.Fields.Type.toString(), b -> b.type = alternateValueSource.Junction.type(),
                NovelSpliceJunctionComparer.Fields.FragmentCount.toString(), b -> b.fragmentCount = alternateValueSource.Junction.fragmentCount(),
                NovelSpliceJunctionComparer.Fields.RegionStart.toString(), b -> b.regionStart = alternateValueSource.Junction.regionStart(),
                NovelSpliceJunctionComparer.Fields.RegionEnd.toString(), b -> b.regionEnd = alternateValueSource.Junction.regionEnd());

        nameToAlternateIndexInitializer = Map.of(
                FLD_GENE_NAME, b -> b.geneName = alternateValueSource.Junction.geneName(),
                FLD_CHROMOSOME, b -> {
                    b.chromosome = alternateValueSource.Junction.chromosome();
                    b.comparisonChromosomeStart = alternateValueSource.ComparisonPositionStart.Chromosome;
                    b.comparisonChromosomeEnd = alternateValueSource.ComparisonPositionEnd.Chromosome;
                },
                FLD_ALT_SJ_POS_START, b -> {
                    b.junctionStart = alternateValueSource.Junction.junctionStart();
                    b.comparisonPositionStart = alternateValueSource.ComparisonPositionStart.Position;
                },
                FLD_ALT_SJ_POS_END, b -> {
                    b.junctionEnd = alternateValueSource.Junction.junctionEnd();
                    b.comparisonPositionEnd = alternateValueSource.ComparisonPositionEnd.Position;
                });
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        // Override since Isofox output is never compared in reportable mode
    }

    @Override
    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Override since Isofox output is never compared in reportable mode
    }

    @Test
    public void fullyMatchesSelfWithLiftoverStart()
    {
        NovelSpliceJunctionData victim = TestNovelSpliceJunctionDataBuilder.BUILDER.create(b -> b.comparisonPositionStart = 5000);
        NovelSpliceJunctionData liftoverVictim = TestNovelSpliceJunctionDataBuilder.BUILDER.create(b -> b.junctionStart = 5000);

        MatchLevel matchLevel = MatchLevel.DETAILED;
        FieldCheckCache fieldConfig = createDefaultThresholds(matchLevel);

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(comparer, liftoverVictim, matchLevel, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(comparer, liftoverVictim, matchLevel, true));
    }

    @Test
    public void fullyMatchesSelfWithLiftoverEnd()
    {
        NovelSpliceJunctionData victim = TestNovelSpliceJunctionDataBuilder.BUILDER.create(b -> b.comparisonPositionEnd = 6000);
        NovelSpliceJunctionData liftoverVictim = TestNovelSpliceJunctionDataBuilder.BUILDER.create(b -> b.junctionEnd = 6000);

        MatchLevel matchLevel = MatchLevel.DETAILED;
        FieldCheckCache fieldConfig = createDefaultThresholds(matchLevel);

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(comparer, liftoverVictim, matchLevel, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(comparer, liftoverVictim, matchLevel, true));
    }

    @Test
    public void keyNonEmptyForLiftoverStart()
    {
        NovelSpliceJunctionData victim = TestNovelSpliceJunctionDataBuilder.BUILDER.create(b -> b.comparisonPositionStart = 5000);
        assertFalse(victim.key().isEmpty());
    }

    @Test
    public void keyNonEmptyForLiftoverEnd()
    {
        NovelSpliceJunctionData victim = TestNovelSpliceJunctionDataBuilder.BUILDER.create(b -> b.comparisonPositionEnd = 6000);
        assertFalse(victim.key().isEmpty());
    }
}
