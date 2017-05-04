package com.hartwig.hmftools.patientreporter.copynumber;

import com.hartwig.hmftools.common.region.GenomeRegion;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Map;

public class CopyNumberAnalysis {

    @NotNull
    private final Map<GenomeRegion, CopyNumberStats> stats;
    @NotNull
    private final List<CopyNumberReport> findings;

    CopyNumberAnalysis(@NotNull final Map<GenomeRegion, CopyNumberStats> stats,
            @NotNull final List<CopyNumberReport> findings) {
        this.stats = stats;
        this.findings = findings;
    }

    @NotNull
    public Map<GenomeRegion, CopyNumberStats> stats() {
        return stats;
    }

    @NotNull
    public List<CopyNumberReport> findings() {
        return findings;
    }
}
