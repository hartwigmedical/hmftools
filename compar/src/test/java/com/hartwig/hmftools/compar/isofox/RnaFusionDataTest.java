package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_DISCORD_FRAGS;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_FILTER;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_KNOWN_TYPE;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_REALIGN_FRAGS;
import static com.hartwig.hmftools.common.rna.RnaFusionFile.FLD_SPLIT_FRAGS;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_JUNC_TYPE_DOWN;
import static com.hartwig.hmftools.compar.isofox.RnaFusionData.FLD_JUNC_TYPE_UP;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchType;

import org.junit.Before;
import org.junit.Test;

public class RnaFusionDataTest extends ComparableItemTest<RnaFusionData, RnaFusionComparer, TestRnaFusionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaFusionComparer(new ComparConfig());
        builder = TestRnaFusionDataBuilder.BUILDER;
        RnaFusionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_KNOWN_TYPE, b -> b.knownType = alternateValueSource.RnaFusion().knownType(),
                FLD_JUNC_TYPE_UP, b -> b.junctionTypeUp = alternateValueSource.RnaFusion().junctionTypeUp(),
                FLD_JUNC_TYPE_DOWN, b -> b.junctionTypeDown = alternateValueSource.RnaFusion().junctionTypeDown(),
                FLD_SPLIT_FRAGS, b -> b.splitFragments = alternateValueSource.RnaFusion().splitFragments(),
                FLD_REALIGN_FRAGS, b -> b.realignedFrags = alternateValueSource.RnaFusion().realignedFrags(),
                FLD_DISCORD_FRAGS, b -> b.discordantFrags = alternateValueSource.RnaFusion().discordantFrags()
        );
        nameToAlternateIndexInitializer = Map.of(
                "FusionName", b -> b.name = alternateValueSource.RnaFusion().name(),
                "ChromosomeUp", b ->
                {
                    b.chromosomeUp = alternateValueSource.RnaFusion().chromosomeUp();
                    b.comparisonChromosomeUp = alternateValueSource.ComparisonPositionUp().Chromosome;
                },
                "ChromosomeDown", b ->
                {
                    b.chromosomeDown = alternateValueSource.RnaFusion().chromosomeDown();
                    b.comparisonChromosomeDown = alternateValueSource.ComparisonPositionDown().Chromosome;
                },
                "PositionUp", b ->
                {
                    b.positionUp = alternateValueSource.RnaFusion().positionUp();
                    b.comparisonPositionUp = alternateValueSource.ComparisonPositionUp().Position;
                },
                "PositionDown", b ->
                {
                    b.positionDown = alternateValueSource.RnaFusion().positionDown();
                    b.comparisonPositionDown = alternateValueSource.ComparisonPositionDown().Position;
                }
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Map.of(
                FLD_FILTER, b -> b.filter = "Filter"
        );
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
    public void fullyMatchesSelfWithLiftoverUp()
    {
        RnaFusionData victim = TestRnaFusionDataBuilder.BUILDER.create(b -> b.comparisonPositionUp = 5000);
        RnaFusionData liftoverVictim = TestRnaFusionDataBuilder.BUILDER.create(b -> b.positionUp = 5000);
        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
    }

    @Test
    public void fullyMatchesSelfWithLiftoverDown()
    {
        RnaFusionData victim = TestRnaFusionDataBuilder.BUILDER.create(b -> b.comparisonPositionDown = 6000);
        RnaFusionData liftoverVictim = TestRnaFusionDataBuilder.BUILDER.create(b -> b.positionDown = 6000);
        DiffThresholds diffThresholds = createDefaultThresholds();

        assertTrue(victim.matches(liftoverVictim));
        assertTrue(liftoverVictim.matches(victim));
        assertNull(victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, false));

        Mismatch expectedMatch = new Mismatch(victim, liftoverVictim, MismatchType.FULL_MATCH, Collections.emptyList());
        assertEquals(expectedMatch, victim.findMismatch(liftoverVictim, MatchLevel.DETAILED, diffThresholds, true));
    }

    @Test
    public void keyNonEmptyForLiftoverUp()
    {
        RnaFusionData victim = TestRnaFusionDataBuilder.BUILDER.create(b -> b.comparisonPositionUp = 5000);
        assertFalse(victim.key().isEmpty());
    }

    @Test
    public void keyNonEmptyForLiftoverDown()
    {
        RnaFusionData victim = TestRnaFusionDataBuilder.BUILDER.create(b -> b.comparisonPositionDown = 6000);
        assertFalse(victim.key().isEmpty());
    }
}