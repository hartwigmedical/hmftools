package com.hartwig.hmftools.patientreporter;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;

import org.jetbrains.annotations.NotNull;

public interface ReportData {

    @NotNull
    List<PatientTumorLocation> patientTumorLocations();

    @NotNull
    Lims limsModel();

    @NotNull
    String signaturePath();

    @NotNull
    String logoRVAPath();

    @NotNull
    String logoCompanyPath();
}
