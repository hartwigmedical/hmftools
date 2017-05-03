package com.hartwig.hmftools.common.copynumber;

import com.hartwig.hmftools.common.slicing.GenomeRegion;
import org.jetbrains.annotations.NotNull;

public class CopyNumber extends GenomeRegion {

    private static final int NORMAL_HUMAN_COPY_NUMBER = 2;

    private final int value;

    public CopyNumber(@NotNull final String chromosome, final long start, final long end, final int value) {
        super(chromosome, start, end);

        assert end >= start;
        assert value >= 0;

        this.value = value;
    }

    public int value() {
        return value;
    }

    public boolean isGain() {
        return value > NORMAL_HUMAN_COPY_NUMBER;
    }

    public boolean isLoss() {
        return value < NORMAL_HUMAN_COPY_NUMBER;
    }
}
