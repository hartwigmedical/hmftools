package com.hartwig.hmftools.patientdb.readers;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.ImmutableTumorTypeLims;
import com.hartwig.hmftools.patientdb.data.TumorTypeLims;

import org.jetbrains.annotations.NotNull;

public class TumorLocationCurationLims {

    @NotNull
    private final Lims lims;
    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    public TumorLocationCurationLims(@NotNull final Lims lims, @NotNull TumorLocationCurator tumorLocationCurator) {
        this.lims = lims;
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    public TumorTypeLims read(@NotNull final String sampleId) {
        final TumorTypeLims limsTumorCurations;

        limsTumorCurations = ImmutableTumorTypeLims.of(sampleId, tumorLocationCurator.search(lims.primaryTumor(sampleId)));

        return limsTumorCurations;
    }

    @NotNull
    public TumorTypeLims readFixedValue(@NotNull final String sampleId) {
        final TumorTypeLims limsTumorCurations;

        limsTumorCurations = ImmutableTumorTypeLims.of(sampleId, tumorLocationCurator.search("Melanoma"));
        return limsTumorCurations;
    }
}
