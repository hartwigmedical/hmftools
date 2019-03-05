package com.hartwig.hmftools.patientdb.readers;

import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.patientdb.LoadClinicalData;
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
    public TumorTypeLims read(@NotNull final String sampleId, @NotNull final String patientId) {
        final TumorTypeLims limsTumorCurations;
        limsTumorCurations = ImmutableTumorTypeLims.of(patientId, tumorLocationCurator.search(lims.primaryTumor(sampleId)));

        return limsTumorCurations;
    }

    @NotNull
    public TumorTypeLims readFixedValue(@NotNull final String patientId) {
        LOGGER.info("COLO829");
        final TumorTypeLims TumorCurationsCOLO;
        LOGGER.info("COLO829");
        TumorCurationsCOLO = ImmutableTumorTypeLims.of(patientId, tumorLocationCurator.search("Melanoma"));
        LOGGER.info("COLO829");
        LOGGER.info(TumorCurationsCOLO);
        return TumorCurationsCOLO;
    }
}
