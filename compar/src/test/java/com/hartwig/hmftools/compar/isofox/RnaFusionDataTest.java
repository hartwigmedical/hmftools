package com.hartwig.hmftools.compar.isofox;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class RnaFusionDataTest extends ComparableItemTest<RnaFusionData, RnaFusionComparer, TestRnaFusionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaFusionComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestRnaFusionDataBuilder.BUILDER;
        RnaFusionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                RnaFusionComparer.Fields.KnownType.toString(), b -> b.knownType = alternateValueSource.KnownType,
                RnaFusionComparer.Fields.JuncTypeUp.toString(), b -> b.junctionTypeUp = alternateValueSource.JunctionTypeUp,
                RnaFusionComparer.Fields.JuncTypeDown.toString(), b -> b.junctionTypeDown = alternateValueSource.JunctionTypeDown,
                RnaFusionComparer.Fields.SplitFrags.toString(), b -> b.splitFragments = alternateValueSource.SplitFragments
        );
        nameToAlternateIndexInitializer = Map.of(
                "FusionName", b -> b.name = alternateValueSource.Name,
                "ChromosomeUp", b ->
                {
                    b.chromosomeUp = alternateValueSource.ChromosomeUp;
                },
                "ChromosomeDown", b ->
                {
                    b.chromosomeDown = alternateValueSource.ChromosomeDown;
                },
                "PositionUp", b ->
                {
                    b.positionUp = alternateValueSource.PositionUp;
                },
                "PositionDown", b ->
                {
                    b.positionDown = alternateValueSource.PositionDown;
                }
        );
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
}
