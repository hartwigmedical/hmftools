package com.hartwig.hmftools.patientdb.readers.drup;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;

import org.jetbrains.annotations.NotNull;

class BaselineReader {

    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    public BaselineReader(@NotNull final TumorLocationCurator tumorLocationCurator) {
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    BaselineData read(@NotNull final EcrfPatient patient) {
        return null;
    }
}
