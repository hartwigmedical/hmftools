package com.hartwig.hmftools.compar.isofox;

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

public class RnaFusionDataTest extends ComparableItemTest<RnaFusionData, RnaFusionComparer, TestRnaFusionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaFusionComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestRnaFusionDataBuilder.BUILDER;
        RnaFusionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                RnaFusionComparer.Fields.KnownFusionType.toString(), b -> b.knownType = alternateValueSource.Fusion.knownType(),
                RnaFusionComparer.Fields.JuncTypeUp.toString(), b -> b.junctionTypeUp = alternateValueSource.Fusion.junctionTypeUp(),
                RnaFusionComparer.Fields.JuncTypeDown.toString(), b -> b.junctionTypeDown = alternateValueSource.Fusion.junctionTypeDown(),
                RnaFusionComparer.Fields.SplitFrags.toString(), b -> b.splitFragments = alternateValueSource.Fusion.splitFragments()
        );
        nameToAlternateIndexInitializer = Map.of(
                "FusionName", b -> b.name = alternateValueSource.Fusion.name(),
                "ChromosomeUp", b ->
                {
                    b.chromosomeUp = alternateValueSource.Fusion.chromosomeUp();
                },
                "ChromosomeDown", b ->
                {
                    b.chromosomeDown = alternateValueSource.Fusion.chromosomeDown();
                },
                "PositionUp", b ->
                {
                    b.positionUp = alternateValueSource.Fusion.positionUp();
                },
                "PositionDown", b ->
                {
                    b.positionDown = alternateValueSource.Fusion.positionDown();
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
