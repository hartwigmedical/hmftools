package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.List;

// Represents a set of possible probes, from which the best one is selected.
public record ProbeCandidateChoice(
        List<ProbeCandidate> candidates,
        ProbeSourceInfo source
)
{
}
