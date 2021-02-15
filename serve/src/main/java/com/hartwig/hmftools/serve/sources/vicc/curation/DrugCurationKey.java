package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;

class DrugCurationKey {

    @NotNull
    private final ViccSource source;
    @NotNull
    private final EvidenceLevel level;
    @NotNull
    private final String treatment;

    public DrugCurationKey(@NotNull final ViccSource source, @NotNull final EvidenceLevel level, @NotNull final String treatment) {
        this.source = source;
        this.level = level;
        this.treatment = treatment;
    }

    @VisibleForTesting
    @NotNull
    ViccSource source() {
        return source;
    }

    @VisibleForTesting
    @NotNull
    EvidenceLevel level() {
        return level;
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
        final DrugCurationKey that = (DrugCurationKey) o;
        return source == that.source && level == that.level && treatment.equals(that.treatment);
    }

    @Override
    public int hashCode() {
        return Objects.hash(source, level, treatment);
    }

    @Override
    public String toString() {
        return "DrugCurationKey{" + "source=" + source + ", evidenceLevel=" + level + ", drugLabels='" + treatment + '\'' + '}';
    }
}
