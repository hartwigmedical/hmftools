package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.List;

// TODO: don't think we need this.
//  The full probe generation and selection should be inside the classes, not done in a central location,
//  because it's generally not possible to express the generation operation as "find all candidates then filter" (e.g. tiling)

// Represents a set of possible probes, from which the best one is selected.
public record ProbeCandidateChoice(
        List<ProbeCandidate> candidates,
        ProbeSourceInfo source
)
{
}
