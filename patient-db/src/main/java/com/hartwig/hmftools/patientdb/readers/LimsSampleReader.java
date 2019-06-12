package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.util.Set;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.data.ImmutableSampleData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LimsSampleReader {

    private static final Logger LOGGER = LogManager.getLogger(LimsSampleReader.class);

    @NotNull
    private final Lims lims;
    @NotNull
    private final Set<String> sequencedSampleIds;

    public LimsSampleReader(@NotNull final Lims lims, @NotNull final Set<String> sequencedSampleIds) {
        this.lims = lims;
        this.sequencedSampleIds = sequencedSampleIds;
    }

    @Nullable
    public SampleData read(@NotNull String sampleId) {
        final LocalDate arrivalDate = lims.arrivalDate(sampleId);
        if (arrivalDate != null) {
            return ImmutableSampleData.of(sampleId,
                    sequencedSampleIds.contains(sampleId),
                    arrivalDate,
                    lims.samplingDate(sampleId),
                    lims.dnaNanograms(sampleId),
                    lims.primaryTumor(sampleId),
                    lims.pathologyTumorPercentage(sampleId));
        } else {
            LOGGER.warn("Skipping sample in LimsSampleReader because arrival date is missing: " + sampleId);
            return null;
        }
    }
}
