package com.hartwig.hmftools.protect;

import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class EventKey {

    @Nullable
    private final String gene;
    @NotNull
    private final String event;

    @NotNull
    public static Set<EventKey> buildUniqueEventSet(@NotNull List<ProtectEvidence> evidences) {
        Set<EventKey> keys = Sets.newHashSet();
        for (ProtectEvidence evidence : evidences) {
            keys.add(create(evidence));
        }
        return keys;
    }

    @NotNull
    public static EventKey create(@NotNull ProtectEvidence evidence) {
        return new EventKey(evidence.gene(), evidence.event());
    }

    private EventKey(@Nullable final String gene, @NotNull final String event) {
        this.gene = gene;
        this.event = event;
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

    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final EventKey that = (EventKey) o;
        return Objects.equals(gene, that.gene) && event.equals(that.event);
    }

    @Override
    public int hashCode() {
        return Objects.hash(gene, event);
    }

    @Override
    public String toString() {
        return gene != null ? gene + " " + event : event;
    }
}
