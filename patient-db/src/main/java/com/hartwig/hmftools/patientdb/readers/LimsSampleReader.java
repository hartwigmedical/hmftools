package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.data.ImmutableSampleData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LimsSampleReader {

    private static final Logger LOGGER = LogManager.getLogger(LimsSampleReader.class);

    @NotNull
    private final Lims lims;

    public LimsSampleReader(@NotNull final Lims lims) {
        this.lims = lims;
    }

    @NotNull
    public List<SampleData> read(@NotNull final Set<String> sampleIds) {
        final List<SampleData> limsBiopsies = Lists.newArrayList();

        sampleIds.forEach(sampleId -> {
            final LocalDate arrivalDate = lims.arrivalDate(sampleId);
            if (arrivalDate != null) {
                limsBiopsies.add(ImmutableSampleData.of(sampleId,
                        arrivalDate,
                        lims.samplingDate(sampleId),
                        lims.dnaNanograms(sampleId),
                        lims.primaryTumor(sampleId),
                        lims.pathologyTumorPercentage(sampleId)));
            } else {
                LOGGER.warn("Skipping sample in LimsSampleReader because arrival date is missing: " + sampleId);
            }
        });
        return limsBiopsies;
    }
}
