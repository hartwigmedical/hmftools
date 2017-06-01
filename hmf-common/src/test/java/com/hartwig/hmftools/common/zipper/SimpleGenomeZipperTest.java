package com.hartwig.hmftools.common.zipper;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Before;
import org.junit.Test;

public class SimpleGenomeZipperTest {

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
        positions = Lists.newArrayList(beforeRegion1, inRegion1, betweenRegion1And2, inRegion2, afterRegion2,
                inRegion3, afterRegion3);
    }

    @Test
    public void testInRegionsOnly() {
        SimpleGenomeZipper.zip(regions, positions, new InRegionHandler());
    }

    @Test
    public void testAllPositions() {
        SimpleGenomeZipper.zip(regions, positions, new AllPositionsHandler());
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

    class AllPositionsHandler implements SimpleGenomeZipperAllPositionsHandler<GenomeRegion, GenomePosition> {

        private int count;

        @Override
        public void handle(@Nullable final GenomeRegion region, @NotNull final GenomePosition position) {
            switch (count) {
                case 0:
                    assertNull(region);
                    assertEquals(beforeRegion1, position);
                    break;
                case 1:
                    assertEquals(region1, region);
                    assertEquals(inRegion1, position);
                    break;
                case 2:
                    assertNull(region);
                    assertEquals(betweenRegion1And2, position);
                    break;
                case 3:
                    assertEquals(region2, region);
                    assertEquals(inRegion2, position);
                    break;
                case 4:
                    assertNull(region);
                    assertEquals(afterRegion2, position);
                    break;
                case 5:
                    assertEquals(region3, region);
                    assertEquals(inRegion3, position);
                    break;
                case 6:
                    assertNull(region);
                    assertEquals(afterRegion3, position);
                    break;
                default:
                    assertTrue(false);
            }
            count++;
        }
    }

    class InRegionHandler implements SimpleGenomeZipperInRegionPositionsHandler<GenomeRegion, GenomePosition> {

        private int count;

        @Override
        public void handle(@NotNull final GenomeRegion region, @NotNull final GenomePosition position) {
            switch (count) {
                case 0:
                    assertEquals(region1, region);
                    assertEquals(inRegion1, position);
                    break;
                case 1:
                    assertEquals(region2, region);
                    assertEquals(inRegion2, position);
                    break;
                case 2:
                    assertEquals(region3, region);
                    assertEquals(inRegion3, position);
                    break;
                default:
                    assertTrue(false);
            }
            count++;
        }
    }
}


