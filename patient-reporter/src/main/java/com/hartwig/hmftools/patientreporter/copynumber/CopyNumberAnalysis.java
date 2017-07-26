package com.hartwig.hmftools.patientreporter.copynumber;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class CopyNumberAnalysis {

    private final int genePanelSize;
    @NotNull
    private final List<CopyNumberReport> findings;

    public CopyNumberAnalysis(final int genePanelSize, @NotNull final List<CopyNumberReport> findings) {
        this.genePanelSize = genePanelSize;
        this.findings = findings;
    }

    public int genePanelSize() {
        return genePanelSize;
    }

    @NotNull
    public List<CopyNumberReport> findings() {
        return findings;
    }
}
