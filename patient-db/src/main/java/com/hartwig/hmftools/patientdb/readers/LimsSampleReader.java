package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class LimsSampleReader {
    private static final Logger LOGGER = LogManager.getLogger(LimsSampleReader.class);

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
                limsBiopsies.add(new SampleData(sampleId, arrivalDate, lims.samplingDateForSample(sampleId)));
            } else {
                LOGGER.warn("Could not find arrival date for sample: " + sampleId);
            }
        });
        return limsBiopsies;
    }
}
