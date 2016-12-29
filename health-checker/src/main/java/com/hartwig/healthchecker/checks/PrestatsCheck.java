package com.hartwig.healthchecker.checks;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

enum PrestatsCheck {
    PRESTATS_PER_BASE_SEQUENCE_QUALITY("Per base sequence quality", true),
    PRESTATS_PER_TILE_SEQUENCE_QUALITY("Per tile sequence quality", true),
    PRESTATS_PER_SEQUENCE_QUALITY_SCORES("Per sequence quality scores", true),
    PRESTATS_PER_BASE_SEQUENCE_CONTENT("Per base sequence content", true),
    PRESTATS_PER_SEQUENCE_GC_CONTENT("Per sequence GC content", true),
    PRESTATS_PER_BASE_N_CONTENT("Per base N content", true),
    PRESTATS_SEQUENCE_LENGTH_DISTRIBUTION("Sequence Length Distribution", true),
    PRESTATS_SEQUENCE_DUPLICATION_LEVELS("Sequence Duplication Levels", true),
    PRESTATS_OVERREPRESENTED_SEQUENCES("Overrepresented sequences", true),
    PRESTATS_ADAPTER_CONTENT("Adapter Content", true),
    PRESTATS_KMER_CONTENT("Kmer Content", true),
    PRESTATS_NUMBER_OF_READS("Total Sequences", false);

    @NotNull
    private final String description;
    private final boolean includeInCountCheck;

    PrestatsCheck(@NotNull final String description, final boolean includeInCountCheck) {
        this.description = description;
        this.includeInCountCheck = includeInCountCheck;
    }

    @NotNull
    public static Optional<PrestatsCheck> getByDescription(@NotNull final String description) {
        return Arrays.stream(PrestatsCheck.values())
                        .filter(prestatsCheck -> prestatsCheck.description.equalsIgnoreCase(description)).findFirst();
    }

    @NotNull
    public static List<PrestatsCheck> valuesToIncludeInCount() {
        List<PrestatsCheck> checksToInclude = Lists.newArrayList();
        for (PrestatsCheck check : values()) {
            if (check.includeInCountCheck) {
                checksToInclude.add(check);
            }
        }
        return checksToInclude;
    }
}
