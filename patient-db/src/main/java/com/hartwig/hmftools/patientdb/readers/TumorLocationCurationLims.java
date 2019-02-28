package com.hartwig.hmftools.patientdb.readers;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.ImmutableTumorTypeLims;
import com.hartwig.hmftools.patientdb.data.TumorTypeLims;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class TumorLocationCurationLims {
    private static final Logger LOGGER = LogManager.getLogger(TumorLocationCurationLims.class);

    @NotNull
    private final Lims lims;
    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    public TumorLocationCurationLims(@NotNull final Lims lims, @NotNull TumorLocationCurator tumorLocationCurator) {
        this.lims = lims;
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    public List<TumorTypeLims> read(@NotNull final List<String> sampleIds) {
        final List<TumorTypeLims> limsTumorCurations = Lists.newArrayList();

        sampleIds.forEach(sampleId -> {
            limsTumorCurations.add(ImmutableTumorTypeLims.of(tumorLocationCurator.search(lims.primaryTumor(sampleId))));

        });
        return limsTumorCurations;
    }

    @NotNull
    public List<TumorTypeLims> readFixedValue(@NotNull final List<String> sampleIds) {
        final List<TumorTypeLims> limsTumorCurations = Lists.newArrayList();

        sampleIds.forEach(sampleId -> {
            limsTumorCurations.add(ImmutableTumorTypeLims.of(tumorLocationCurator.search("Melanoma")));
        });
        return limsTumorCurations;
    }
}
