package com.hartwig.hmftools.common.protect;

import java.util.Comparator;
import java.util.Objects;

import org.apache.commons.lang3.StringUtils;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EvidenceComparator implements Comparator<ProtectEvidence> {

    @Override
    public int compare(@NotNull ProtectEvidence evidence1, @NotNull ProtectEvidence evidence2) {
        int reportedCompare = -Boolean.compare(evidence1.reported(), evidence2.reported());
        if (reportedCompare != 0) {
            return reportedCompare;
        }

        int geneCompare = StringUtils.compare(evidence1.gene(), evidence2.gene());
        if (geneCompare != 0) {
            return geneCompare;
        }

        int transcriptCompare = StringUtils.compare(evidence1.transcript(), evidence2.transcript());
        if (transcriptCompare != 0) {
            return transcriptCompare;
        }

        int isCanonicalCompare = compareBoolean(evidence1.isCanonical(), evidence2.isCanonical());
        if (isCanonicalCompare != 0) {
            return isCanonicalCompare;
        }

        int eventCompare = evidence1.event().compareTo(evidence2.event());
        if (eventCompare != 0) {
            return eventCompare;
        }

        int levelCompare = evidence1.level().compareTo(evidence2.level());
        if (levelCompare != 0) {
            return levelCompare;
        }

        int onLabelCompare = -Boolean.compare(evidence1.onLabel(), evidence2.onLabel());
        if (onLabelCompare != 0) {
            return onLabelCompare;
        }

        int treatmentCompare = evidence1.treatment().compareTo(evidence2.treatment());
        if (treatmentCompare != 0) {
            return treatmentCompare;
        }

        int directionCompare = evidence1.direction().compareTo(evidence2.direction());
        if (directionCompare != 0) {
            return directionCompare;
        }

        return 0;
    }

    private static int compareBoolean(@Nullable Boolean boolean1, @Nullable Boolean boolean2) {
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
