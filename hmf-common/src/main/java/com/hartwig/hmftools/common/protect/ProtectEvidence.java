package com.hartwig.hmftools.common.protect;

import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.commons.lang3.StringUtils;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ProtectEvidence implements Comparable<ProtectEvidence> {

    @Nullable
    public abstract String gene();

    @Nullable
    public abstract String transcript();

    @Nullable
    public abstract Boolean isCanonical();

    @NotNull
    public abstract String event();

    @Nullable
    public abstract Boolean eventIsHighDriver();

    public abstract boolean germline();

    public abstract boolean reported();

    @NotNull
    public abstract String treatment();

    public abstract boolean onLabel();

    @NotNull
    public abstract EvidenceLevel level();

    @NotNull
    public abstract EvidenceDirection direction();

    @NotNull
    public abstract Set<ProtectSource> protectSources();

    @Override
    public int compareTo(@NotNull final ProtectEvidence o) {
        int reportedCompare = -Boolean.compare(reported(), o.reported());
        if (reportedCompare != 0) {
            return reportedCompare;
        }

        int geneCompare = StringUtils.compare(gene(), o.gene());
        if (geneCompare != 0) {
            return geneCompare;
        }

        int transcriptCompare = StringUtils.compare(transcript(), o.transcript());
        if (transcriptCompare != 0) {
            return transcriptCompare;
        }

        int isCanonicalCompare = compareBoolean(isCanonical(), o.isCanonical());
        if (isCanonicalCompare != 0) {
            return isCanonicalCompare;
        }

        int eventCompare = event().compareTo(o.event());
        if (eventCompare != 0) {
            return eventCompare;
        }

        int levelCompare = level().compareTo(o.level());
        if (levelCompare != 0) {
            return levelCompare;
        }

        int onLabelCompare = -Boolean.compare(onLabel(), o.onLabel());
        if (onLabelCompare != 0) {
            return onLabelCompare;
        }

        int treatmentCompare = treatment().compareTo(o.treatment());
        if (treatmentCompare != 0) {
            return treatmentCompare;
        }

        int directionCompare = direction().compareTo(o.direction());
        if (directionCompare != 0) {
            return directionCompare;
        }

        return 0;
    }

    public static int compareBoolean(@Nullable Boolean boolean1, @Nullable Boolean boolean2) {
        if (Objects.equals(boolean1, boolean2)) {
            return 0;
        } else if (boolean1 == null) {
            return -1;
        } else if (boolean2 == null) {
            return 1;
        } else {
            return boolean1.compareTo(boolean2);
        }
    }
}