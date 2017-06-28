package com.hartwig.hmftools.common.zipper;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import mockit.Expectations;
import mockit.Injectable;

public class SimpleGenomeZipperTest {

    @Injectable
    private SimpleGenomeZipperAllPositionsHandler<GenomeRegion, GenomePosition> allPositionsHandler;

    @Injectable
    private SimpleGenomeZipperInRegionPositionsHandler<GenomeRegion, GenomePosition> inRegionPositionsHandler;

    private List<GenomeRegion> regions;
    private List<GenomePosition> positions;

    private GenomeRegion region1;
    private GenomeRegion region2;
    private GenomeRegion region3;

    private GenomePosition beforeRegion1;
    private GenomePosition inRegion1;
    private GenomePosition betweenRegion1And2;
    private GenomePosition inRegion2;
    private GenomePosition afterRegion2;
    private GenomePosition inRegion3;
    private GenomePosition afterRegion3;

    @Before
    public void setup() {

        region1 = createRegion("1", 200, 300);
        region2 = createRegion("1", 400, 500);
        region3 = createRegion("2", 1, 100);
        regions = Lists.newArrayList(region1, region2, region3);

        beforeRegion1 = createPosition("1", 50);
        inRegion1 = createPosition("1", 250);
        betweenRegion1And2 = createPosition("1", 350);
        inRegion2 = createPosition("1", 450);
        afterRegion2 = createPosition("1", 550);
        inRegion3 = createPosition("2", 50);
        afterRegion3 = createPosition("2", 150);
        positions = Lists.newArrayList(beforeRegion1, inRegion1, betweenRegion1And2, inRegion2, afterRegion2, inRegion3,
                afterRegion3);
    }

    @Test
    public void testInRegionsOnly() {
        new Expectations() {{
            inRegionPositionsHandler.handle(region1, inRegion1);
            inRegionPositionsHandler.handle(region2, inRegion2);
            inRegionPositionsHandler.handle(region3, inRegion3);
        }};

        SimpleGenomeZipper.zip(regions, positions, inRegionPositionsHandler);
    }

    @Test
    public void testAllPositions() {
        new Expectations() {{
            allPositionsHandler.handle(null, beforeRegion1);
            allPositionsHandler.handle(region1, inRegion1);
            allPositionsHandler.handle(null, betweenRegion1And2);
            allPositionsHandler.handle(region2, inRegion2);
            allPositionsHandler.handle(null, afterRegion2);
            allPositionsHandler.handle(region3, inRegion3);
            allPositionsHandler.handle(null, afterRegion3);
        }};

        SimpleGenomeZipper.zip(regions, positions, allPositionsHandler);
    }

    private GenomeRegion createRegion(final String chromosome, final long start, final long end) {
        return ImmutableBEDGenomeRegion.builder().chromosome(chromosome).start(start).end(end).build();
    }

    private GenomePosition createPosition(final String chromosome, final long position) {
        return new GenomePosition() {
            @NotNull
            @Override
            public String chromosome() {
                return chromosome;
            }

            @Override
            public long position() {
                return position;
            }
        };
    }
}


