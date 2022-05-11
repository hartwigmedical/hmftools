package com.hartwig.hmftools.protect;

import java.util.Objects;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class EvidenceKey {

    @Nullable
    private final String gene;
    @NotNull
    private final String event;
    @NotNull
    private final String treatment;

    @NotNull
    public static Set<EvidenceKey> buildUniqueEventSet(@NotNull Iterable<ProtectEvidence> evidences) {
        Set<EvidenceKey> keys = Sets.newHashSet();
        for (ProtectEvidence evidence : evidences) {
            keys.add(create(evidence));
        }
        return keys;
    }

    @NotNull
    public static EvidenceKey create(@NotNull ProtectEvidence evidence) {
        return new EvidenceKey(evidence.gene(), evidence.event(), evidence.treatment());
    }

    private EvidenceKey(@Nullable final String gene, @NotNull final String event, @NotNull final String treatment) {
        this.gene = gene;
        this.event = event;
        this.treatment = treatment;
    }

    @VisibleForTesting
    @Nullable
    String gene() {
        return gene;
    }

    @VisibleForTesting
    @NotNull
    String event() {
        return event;
    }

    @VisibleForTesting
    @NotNull
    String treatment() {
        return treatment;
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final EvidenceKey evidenceKey = (EvidenceKey) o;
        return Objects.equals(gene, evidenceKey.gene) && event.equals(evidenceKey.event) && treatment.equals(evidenceKey.treatment);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, event, treatment);
    }

    @Override
    public String toString() {
        return "EventKey{" + "gene='" + gene + '\'' + ", event='" + event + '\'' + ", treatment='" + treatment + '\'' + '}';
    }
}
