package com.hartwig.hmftools.svanalysis.svgraph;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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
}
