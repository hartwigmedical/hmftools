package com.hartwig.healthchecker.common.flagstatreader;

import java.util.Arrays;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public enum FlagStatsType {
    TOTAL_INDEX(0),
    SECONDARY_INDEX(1),
    SUPPLEMENTARY_INDEX(2),
    DUPLICATES_INDEX(3),
    MAPPED_INDEX(4),
    PAIRED_IN_SEQ_INDEX(5),
    READ_1_INDEX(6),
    READ_2_INDEX(7),
    PROPERLY_PAIRED_INDEX(8),
    ITSELF_AND_MATE_INDEX(9),
    SINGLETONS_INDEX(10),
    MATE_MAP_DIF_CHR_INDEX(11),
    MATE_MAP_DIF_CHR_Q5_INDEX(12);

    private final int index;

    FlagStatsType(final int index) {
        this.index = index;
    }

    @NotNull
    public static Optional<FlagStatsType> getByIndex(final int index) {
        return Arrays.stream(FlagStatsType.values()).filter(type -> type.getIndex() == index).findFirst();
    }

    public int getIndex() {
        return index;
    }
}
