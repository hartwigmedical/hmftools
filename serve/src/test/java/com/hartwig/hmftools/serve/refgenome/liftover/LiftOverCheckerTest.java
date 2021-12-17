package com.hartwig.hmftools.serve.refgenome.liftover;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LiftOverCheckerTest {

    @Test
    public void canCheckLiftOverPositions() {
        GenomePosition position1 = new TestGenomePosition("1", 1000);
        LiftOverResult lifted = ImmutableLiftOverResult.builder().chromosome("1").position(1000).build();
        assertTrue(LiftOverChecker.isValidLiftedPosition(lifted, position1));
        assertFalse(LiftOverChecker.isValidLiftedPosition(null, position1));

        GenomePosition position2 = new TestGenomePosition("X", 1314965);
        assertFalse(LiftOverChecker.isValidLiftedPosition(null, position2));
    }

    @Test
    public void canCheckLiftOverRegions() {
        GenomeRegion region1 = new TestGenomeRegion("1", 1000, 2000);
        LiftOverResult liftedStart = ImmutableLiftOverResult.builder().chromosome("1").position(1000).build();
        LiftOverResult liftedEnd = ImmutableLiftOverResult.builder().chromosome("1").position(2000).build();
        assertTrue(LiftOverChecker.isValidLiftedRegion(liftedStart, liftedEnd, region1));
        assertFalse(LiftOverChecker.isValidLiftedRegion(null, liftedEnd, region1));
        assertFalse(LiftOverChecker.isValidLiftedRegion(liftedStart, null, region1));

        GenomeRegion region2 = new TestGenomeRegion("X", 1314965, 1314966);
        assertFalse(LiftOverChecker.isValidLiftedRegion(null, null, region2));
    }

    private static class TestGenomePosition implements GenomePosition {

        @NotNull
        private final String chromosome;
        private final int position;

        public TestGenomePosition(@NotNull final String chromosome, final int position) {
            this.chromosome = chromosome;
            this.position = position;
        }

        @Override
        @NotNull
        public String chromosome() {
            return chromosome;
        }

        @Override
        public int position() {
            return position;
        }

        @Override
        public String toString() {
            return "GenomePositionImpl{" + "chromosome='" + chromosome + '\'' + ", position=" + position + '}';
        }
    }

    private static class TestGenomeRegion implements GenomeRegion {

        @NotNull
        private final String chromosome;
        private final int start;
        private final int end;

        public TestGenomeRegion(@NotNull final String chromosome, final int start, final int end) {
            this.chromosome = chromosome;
            this.start = start;
            this.end = end;
        }

        @NotNull
        @Override
        public String chromosome() {
            return chromosome;
        }

        @Override
        public int start() {
            return start;
        }

        @Override
        public int end() {
            return end;
        }
    }
}