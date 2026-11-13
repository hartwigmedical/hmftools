package com.hartwig.hmftools.compar.isofox;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class RnaFusionDataTest extends ComparableItemTest<RnaFusionData, RnaFusionComparer, TestRnaFusionDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaFusionComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestRnaFusionDataBuilder.BUILDER;
        RnaFusionData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                RnaFusionComparer.Fields.JuncTypeUp.toString(), b -> b.junctionTypeUp = alternateValueSource.JunctionTypeUp,
                RnaFusionComparer.Fields.JuncTypeDown.toString(), b -> b.junctionTypeDown = alternateValueSource.JunctionTypeDown,
                RnaFusionComparer.Fields.SplitFrags.toString(), b -> b.splitFragments = alternateValueSource.SplitFragments,
                RnaFusionComparer.Fields.Filter.toString(), b -> b.filter = alternateValueSource.Filter
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
        reportabilityFieldToFalseReportabilityInitializer = Map.of(
                RnaFusionComparer.Fields.KnownType.toString(), b -> b.knownType = KnownFusionType.NONE);
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
