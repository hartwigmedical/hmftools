package com.hartwig.hmftools.compar.isofox;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.KnownFusionType;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestRnaFusionDataBuilder
{
    public String name = "TMPRSS2_ERG";
    public String chromosomeUp = "21";
    public int positionUp = 1000;
    public String chromosomeDown = "21";
    public int positionDown = 2000;
    public KnownFusionType knownType = KnownFusionType.KNOWN_PAIR;
    public String junctionTypeUp = "KNOWN";
    public String junctionTypeDown = "KNOWN";
    public int splitFragments = 100;
    public int realignedFrags = 50;
    public int discordantFrags = 20;
    public String comparisonChromosomeUp = "21";
    public int comparisonPositionUp = 1000;
    public String comparisonChromosomeDown = "21";
    public int comparisonPositionDown = 2000;
    public String filter = "PASS";

    private static final Consumer<TestRnaFusionDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.name = "EML4_ALK";
        b.chromosomeUp = "2";
        b.positionUp = 3000;
        b.chromosomeDown = "2";
        b.positionDown = 4000;
        b.knownType = KnownFusionType.KNOWN_PROM3;
        b.junctionTypeUp = "CANONICAL";
        b.junctionTypeDown = "CANONICAL";
        b.splitFragments = 10;
        b.realignedFrags = 5;
        b.discordantFrags = 2;
        b.comparisonChromosomeUp = "2";
        b.comparisonPositionUp = 3000;
        b.comparisonChromosomeDown = "2";
        b.comparisonPositionDown = 4000;
        b.filter = "PASS";
    };

    public static final TestComparableItemBuilder<TestRnaFusionDataBuilder, RnaFusionData> BUILDER =
            new TestComparableItemBuilder<>(TestRnaFusionDataBuilder::new, TestRnaFusionDataBuilder::build, ALTERNATE_INITIALIZER);

    private RnaFusionData build()
    {
        final RnaFusion fusion = ImmutableRnaFusion.builder()
                .name(name)
                .chromosomeUp(chromosomeUp)
                .chromosomeDown(chromosomeDown)
                .positionUp(positionUp)
                .positionDown(positionDown)
                .orientationUp((byte) 1)
                .orientationDown((byte) -1)
                .junctionTypeUp(junctionTypeUp)
                .junctionTypeDown(junctionTypeDown)
                .knownType(knownType)
                .svType(StructuralVariantType.BND)
                .splitFragments(splitFragments)
                .realignedFrags(realignedFrags)
                .discordantFrags(discordantFrags)
                .depthUp(-1)
                .depthDown(-1)
                .maxAnchorLengthUp(-1)
                .maxAnchorLengthDown(-1)
                .cohortFrequency(-1)
                .filter(filter)
                .build();

        return new RnaFusionData(fusion,
                new BasePosition(comparisonChromosomeUp, comparisonPositionUp),
                new BasePosition(comparisonChromosomeDown, comparisonPositionDown));
    }
}