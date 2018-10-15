package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.List;

public class SimpleSimplificationStrategy implements SimplificationStrategy {
    private static final Logger LOGGER = LogManager.getLogger(BreakpointGraph.class);
    private static final double MAX_COPY_NUMBER_INCONSISTENCY = 0.5;
    private static final double MAX_SV_PLOIDY_INCONSISTENCY = 0.75;
    @Override
    public boolean shouldSimplify(Simplification s) {
        for (BreakendConsistency bc : s.consistency()) {
            switch (s.type()) {
                case SimpleDeletion:
                case SimpleDuplication:
                    if (bc.copyNumberDelta() > MAX_COPY_NUMBER_INCONSISTENCY) {
                        LOGGER.debug("Not simplifying {}: ∆cn={} ploidy={}", s.type(), bc.copyNumberDelta(), bc.ploidy());
                        return false;
                    }
                    break;
                case Chain:
                    break;
            }
            switch (s.type()) {
                case SimpleDeletion:
                case SimpleDuplication:
                case Chain:
                    if (bc.eventDelta() > MAX_SV_PLOIDY_INCONSISTENCY) {
                        LOGGER.debug("Not simplifying {}: ∆sv={} ploidy={}", s.type(), bc.eventDelta(), bc.ploidy());
                        return false;
                    }
                    break;
            }
        }
        return true;
    }

    @Override
    public boolean couldBeDirectlyLinked(EnrichedStructuralVariant sv, BgSegment segment) {
        if (sv.ploidy() < segment.maxMajorPloidy() + MAX_COPY_NUMBER_INCONSISTENCY) {
            return true;
        } else {
            // we can't traverse across this segment since our CN is too low
            LOGGER.debug("Unable to traverse across {} from {}. ploidy={}, maxMajorCN={}", segment, sv, sv.ploidy(), segment.maxMajorPloidy());
            return false;
        }
    }

    @Override
    public boolean couldBeDirectlyLinked(
            EnrichedStructuralVariant sv1,
            EnrichedStructuralVariant sv2) {
        double diff = sv1.ploidy() - sv2.ploidy();
        if (Math.abs(diff) < MAX_SV_PLOIDY_INCONSISTENCY) {
            return true;
        } else {
            // we can't traverse across this segment since our CN is too low
            LOGGER.debug("Unable to traverse from {} to {}. ∆ploidy={}", sv1, sv2, diff);
            return false;
        }
    }

    @Override
    public boolean couldBeFoldBackLinked(EnrichedStructuralVariant foldback, List<BgSegment> segments) {
        return segments.stream().allMatch(segment -> foldback.ploidy() < segment.maxMajorPloidy() + MAX_COPY_NUMBER_INCONSISTENCY);
    }

    @Override
    public boolean couldBeFoldBackLinked(
            EnrichedStructuralVariant sv1,
            EnrichedStructuralVariant sv2) {
        if (sv1.isFoldBackInversion()) {
            double diff = sv1.ploidy() - sv2.ploidy() / 2;
            if (Math.abs(diff) < MAX_SV_PLOIDY_INCONSISTENCY) {
                return true;
            }
        }
        if (sv2.isFoldBackInversion()) {
            double diff = sv1.ploidy() / 2 - sv2.ploidy();
            if (Math.abs(diff) < MAX_SV_PLOIDY_INCONSISTENCY) {
                return true;
            }
        }
        return false;
    }
}
