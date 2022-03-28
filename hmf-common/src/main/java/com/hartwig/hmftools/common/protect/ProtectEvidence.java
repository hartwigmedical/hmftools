package com.hartwig.hmftools.common.protect;

import java.util.Objects;
import java.util.Set;

import com.hartwig.hmftools.common.serve.Knowledgebase;
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

    @NotNull
    public abstract String event();

    @Nullable
    public abstract Boolean eventIsHighDriver();

    @NotNull
    public abstract ProtectEvidenceType evidenceType();

    @Nullable
    public abstract Integer rangeRank();

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
    public abstract Set<String> evidenceUrls();

    @NotNull
    public abstract Set<Knowledgebase> sources();

    @NotNull
    public abstract String sourceEvent();

    @NotNull
    public abstract Set<String> sourceUrls();


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

        int eventCompare = event().compareTo(o.event());
        if (eventCompare != 0) {
            return eventCompare;
        }

        int evidenceTypeCompare = evidenceType().compareTo(o.evidenceType());
        if (evidenceTypeCompare != 0) {
            return evidenceTypeCompare;
        }

        int rangeRankCompare = compareInteger(rangeRank(), o.rangeRank());
        if (rangeRankCompare != 0) {
            return rangeRankCompare;
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

        int sourceEventCompare = sourceEvent().compareTo(o.sourceEvent());
        if (sourceEventCompare != 0) {
            return sourceEventCompare;
        }

        int germlineCompare = -Boolean.compare(germline(), o.germline());
        if (germlineCompare != 0) {
            return germlineCompare;
        }

        int isHighDriverCompare = compareBoolean(eventIsHighDriver(), o.eventIsHighDriver());
        if (isHighDriverCompare != 0) {
            return isHighDriverCompare;
        }
        return 0;
    }

    public static int compareInteger(@Nullable Integer int1, @Nullable Integer int2) {
        if (Objects.equals(int1, int2)) {
            return 0;
        } else if (int1 == null) {
            return -1;
        } else if (int2 == null) {
            return 1;
        } else {
            return int1.compareTo(int2);
        }
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
