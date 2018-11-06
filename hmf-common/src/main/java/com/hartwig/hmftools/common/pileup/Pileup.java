package com.hartwig.hmftools.common.pileup;

import java.util.Map;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Pileup extends GenomePosition {

    int readCount();

    @NotNull
    String referenceBase();

    int referenceCount();

    int gMismatchCount();

    int aMismatchCount();

    int tMismatchCount();

    int cMismatchCount();

    Map<String, Integer> insertionCounts();

    Map<String, Integer> deletionCounts();

    default int insertions() {
        return insertionCounts().values().stream().mapToInt(x -> x).sum();
    }

    default int inframeInsertions() {
        return insertionCounts().entrySet().stream().filter(x -> (x.getKey().length() - 1) % 3 == 0).mapToInt(Map.Entry::getValue).sum();
    }

    default int deletions() {
        return deletionCounts().values().stream().mapToInt(x -> x).sum();
    }

    default int inframeDeletions() {
        return deletionCounts().entrySet().stream().filter(x -> (x.getKey().length() - 1) % 3 == 0).mapToInt(Map.Entry::getValue).sum();
    }

    default int indels() {
        return insertions() + deletions();
    }

    default int inframeIndels() {
        return inframeInsertions() + inframeDeletions();
    }

    default int mismatchCount(char base) {
        switch (Character.toUpperCase(base)) {
            case 'G':
                return gMismatchCount();
            case 'A':
                return aMismatchCount();
            case 'T':
                return tMismatchCount();
            case 'C':
                return cMismatchCount();
            default:
                return 0;
        }
    }
}
