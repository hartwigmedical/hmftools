package com.hartwig.hmftools.patientdb.readers.wide;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.curators.BiopsySiteCurator;
import com.hartwig.hmftools.patientdb.data.BiopsyData;

import org.jetbrains.annotations.NotNull;

public class BiopsyReader {

    @NotNull
    private final BiopsySiteCurator biopsySiteCurator;

    BiopsyReader(@NotNull final BiopsySiteCurator biopsySiteCurator) {
        this.biopsySiteCurator = biopsySiteCurator;
    }

    @NotNull
    List<BiopsyData> read() {
        final List<BiopsyData> biopsies = Lists.newArrayList();
        return biopsies;
    }
}
