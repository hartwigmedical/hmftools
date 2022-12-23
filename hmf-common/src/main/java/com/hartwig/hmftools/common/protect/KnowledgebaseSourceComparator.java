package com.hartwig.hmftools.common.protect;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.Comparator;
import java.util.Objects;

public class KnowledgebaseSourceComparator implements Comparator<KnowledgebaseSource> {

    @Override
    public int compare(@NotNull KnowledgebaseSource source1, @NotNull KnowledgebaseSource source2) {

        int knowledgebaseNameCompare = source1.name().toString().compareTo(source2.name().toString());
        if (knowledgebaseNameCompare != 0) {
            return knowledgebaseNameCompare;
        }

        int sourceEventCompare = source1.sourceEvent().compareTo(source2.sourceEvent());
        if (sourceEventCompare != 0) {
            return sourceEventCompare;
        }

        int evidenceTypeCompare = source1.evidenceType().display().compareTo(source2.evidenceType().display());
        if (evidenceTypeCompare != 0) {
            return evidenceTypeCompare;
        }

        int rangerankCompare = compareInteger(source1.rangeRank(), source2.rangeRank());
        if (rangerankCompare != 0) {
            return rangerankCompare;
        }

        return 0;
    }

    private static int compareInteger(@Nullable Integer integer1, @Nullable Integer integer2) {
        if (Objects.equals(integer1, integer2)) {
            return 0;
        } else if (integer1 == null) {
            return -1;
        } else if (integer2 == null) {
            return 1;
        } else {
            return integer1.compareTo(integer2);
        }
    }
}