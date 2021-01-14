package com.hartwig.hmftools.protect.evidence;

import java.util.Set;

import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public class PersonalizedEvidenceFactory {

    @NotNull
    private final Set<String> patientTumorDoids;

    public PersonalizedEvidenceFactory(@NotNull final Set<String> patientTumorDoids) {
        this.patientTumorDoids = patientTumorDoids;
    }

    @NotNull
    public ImmutableProtectEvidence.Builder somaticallyReportableEvidence(@NotNull ActionableEvent event) {
        return patientEvidenceBuilder(event).reported(true).germline(false);
    }

    @NotNull
    public ImmutableProtectEvidence.Builder patientEvidenceBuilder(@NotNull ActionableEvent actionable) {
        return ImmutableProtectEvidence.builder()
                .treatment(actionable.treatment())
                .level(actionable.level())
                .direction(actionable.direction())
                .onLabel(patientTumorDoids.contains(actionable.doid()))
                .addSources(actionable.source())
                .urls(actionable.urls());
    }
}
