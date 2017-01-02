package com.hartwig.hmftools.patientreporter;

import java.util.List;

import com.google.common.base.Predicates;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.jetbrains.annotations.NotNull;

public class ConsensusRule {

    @NotNull
    private final Slicer highConfidenceSlicer;
    @NotNull
    private final Slicer cpctGenePanelSlicer;

    public ConsensusRule(@NotNull final Slicer highConfidenceSlicer, @NotNull final Slicer cpctGenePanelSlicer) {
        this.highConfidenceSlicer = highConfidenceSlicer;
        this.cpctGenePanelSlicer = cpctGenePanelSlicer;
    }

    @NotNull
    public List<SomaticVariant> applyConsensusRule(@NotNull List<SomaticVariant> variants) {
        Predicates.and();
        return null;

    }
}
