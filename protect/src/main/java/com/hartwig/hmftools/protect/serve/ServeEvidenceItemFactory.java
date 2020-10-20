package com.hartwig.hmftools.protect.serve;

import java.util.Set;

import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public class ServeEvidenceItemFactory {

    @NotNull
    public static ServeEvidenceItem create(@NotNull String genomicEvent, @NotNull Set<String> doid, @NotNull ActionableEvent actionable) {
        return builder(doid, actionable).genomicEvent(genomicEvent).reported(true).build();
    }

    @NotNull
    public static ImmutableServeEvidenceItem.Builder builder(@NotNull Set<String> doid, @NotNull ActionableEvent actionable) {
        return ImmutableServeEvidenceItem.builder()
                .source(actionable.source())
                .direction(actionable.direction())
                .level(actionable.level())
                .treatment(actionable.treatment())
                .onLabel(doid.contains(actionable.doid()))
                .direction(actionable.direction())
                .url(actionable.url());
    }

}
