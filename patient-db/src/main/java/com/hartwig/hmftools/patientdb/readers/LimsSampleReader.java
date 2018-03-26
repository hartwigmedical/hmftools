package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.data.ImmutableSampleData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.jetbrains.annotations.NotNull;

class LimsSampleReader {
    @NotNull
    private final Lims lims;

    LimsSampleReader(@NotNull final Lims lims) {
        this.lims = lims;
    }

    @NotNull
    List<SampleData> read(@NotNull final List<String> sampleIds) {
        final List<SampleData> limsBiopsies = Lists.newArrayList();

        sampleIds.forEach(sampleId -> {
            final LocalDate arrivalDate = lims.arrivalDateForSample(sampleId);
            if (arrivalDate != null) {
                limsBiopsies.add(ImmutableSampleData.of(sampleId, arrivalDate, lims.samplingDateForSample(sampleId),
                        lims.dnaNanogramsForSample(sampleId),
                        lims.tumorPercentageForSample(sampleId)));
            }
        });
        return limsBiopsies;
    }
}
