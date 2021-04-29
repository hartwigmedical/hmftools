package com.hartwig.hmftools.serve.refgenome.liftover;

import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class LiftOverChecker {

    private static final Logger LOGGER = LogManager.getLogger(LiftOverChecker.class);
    private static final Set<GenomePosition> UNLIFTABLE_POSITIONS_37_TO_38 = Sets.newHashSet();

    static {
        UNLIFTABLE_POSITIONS_37_TO_38.add(new GenomePositionImpl("X", 1314965));
        UNLIFTABLE_POSITIONS_37_TO_38.add(new GenomePositionImpl("X", 1314966));
    }

    private LiftOverChecker() {
    }

    public static boolean isValidLiftedPosition(@Nullable LiftOverResult liftedResult, @NotNull GenomePosition liftedPosition) {
        if (liftedResult == null) {
            GenomePosition converted = new GenomePositionImpl(liftedPosition.chromosome(), liftedPosition.position());
            String message = String.format(" Liftover could not be performed on '%s'", liftedPosition);
            if (UNLIFTABLE_POSITIONS_37_TO_38.contains(converted)) {
                LOGGER.debug(message);
            } else {
                LOGGER.warn(message);
            }
            return false;
        }

        return true;
    }

    public static boolean isValidLiftedRegion(@Nullable LiftOverResult liftedStart, @Nullable LiftOverResult liftedEnd,
            @NotNull GenomeRegion liftedRegion) {
        return isValidLiftedPosition(liftedStart, new GenomePositionImpl(liftedRegion.chromosome(), liftedRegion.start())) &&
                isValidLiftedPosition(liftedEnd, new GenomePositionImpl(liftedRegion.chromosome(), liftedRegion.end()));
    }

    private static class GenomePositionImpl implements GenomePosition {

        @NotNull
        private final String chromosome;
        private final long position;

        public GenomePositionImpl(@NotNull final String chromosome, final long position) {
            this.chromosome = chromosome;
            this.position = position;
        }

        @Override
        @NotNull
        public String chromosome() {
            return chromosome;
        }

        @Override
        public long position() {
            return position;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final GenomePositionImpl that = (GenomePositionImpl) o;
            return position == that.position && chromosome.equals(that.chromosome);
        }

        @Override
        public int hashCode() {
            return Objects.hash(chromosome, position);
        }

        @Override
        public String toString() {
            return "GenomePositionImpl{" + "chromosome='" + chromosome + '\'' + ", position=" + position + '}';
        }
    }
}
