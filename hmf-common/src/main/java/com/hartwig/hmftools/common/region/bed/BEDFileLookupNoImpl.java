package com.hartwig.hmftools.common.region.bed;

import java.io.IOException;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class BEDFileLookupNoImpl implements BEDFileLookup {
    @Override
    public double score(@NotNull final GenomePosition position) throws IOException {
        return 0;
    }

    @Override
    public void close() throws IOException {

    }
}
