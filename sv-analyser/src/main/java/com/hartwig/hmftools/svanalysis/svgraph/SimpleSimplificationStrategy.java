package com.hartwig.hmftools.svanalysis.svgraph;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SimpleSimplificationStrategy implements SimplificationStrategy {
    private static final Logger LOGGER = LogManager.getLogger(BreakpointGraph.class);
    private static final double MAX_COPY_NUMBER_INCONSISTENCY = 0.5;
    private static final double MAX_SV_PLOIDY_INCONSISTENCY = 0.5;
    @Override
    public boolean shouldSimplify(Simplification s) {
        for (BreakendConsistency bc : s.consistency()) {
            if (bc.copyNumberDelta() > MAX_COPY_NUMBER_INCONSISTENCY) {
                LOGGER.debug("Not simplifying {}: ∆cn={}", s.type(), bc.copyNumberDelta());
                return false;
            }
            if (bc.eventDelta() > MAX_SV_PLOIDY_INCONSISTENCY) {
                LOGGER.debug("Not simplifying {}: ∆sv={}", s.type(), bc.copyNumberDelta());
                return false;
            }
        }
        return true;
    }
}
