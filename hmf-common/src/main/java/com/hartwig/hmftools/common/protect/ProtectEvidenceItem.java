package com.hartwig.hmftools.common.protect;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ProtectEvidenceItem implements Comparable<ProtectEvidenceItem> {

    @NotNull
    public abstract String genomicEvent();

    @NotNull
    public abstract Knowledgebase source();

    public abstract boolean reported();

    @NotNull
    public abstract String treatment();

    public abstract boolean onLabel();

    @NotNull
    public abstract EvidenceLevel level();

    @NotNull
    public abstract EvidenceDirection direction();

    @NotNull
    public abstract String url();

    @Override
    public int compareTo(@NotNull final ProtectEvidenceItem o) {

        int reportedCompare = -Boolean.compare(reported(), o.reported());
        if (reportedCompare != 0) {
            return reportedCompare;
        }

        int eventCompare = genomicEvent().compareTo(o.genomicEvent());
        if (eventCompare != 0) {
            return eventCompare;
        }

        int treatmentCompare = treatment().compareTo(o.treatment());
        if (treatmentCompare != 0) {
            return treatmentCompare;
        }

        int onLabelCompare = -Boolean.compare(onLabel(), o.onLabel());
        if (onLabelCompare != 0) {
            return onLabelCompare;
        }

        int levelCompare = level().compareTo(o.level());
        if (levelCompare != 0) {
            return levelCompare;
        }

        int directionCompare = direction().compareTo(o.direction());
        if (directionCompare != 0) {
            return directionCompare;
        }

        int sourceCompare = source().compareTo(o.source());
        if (sourceCompare != 0) {
            return directionCompare;
        }

        return url().compareTo(o.url());
    }
}
