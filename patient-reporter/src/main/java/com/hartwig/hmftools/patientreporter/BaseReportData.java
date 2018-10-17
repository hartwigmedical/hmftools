package com.hartwig.hmftools.patientreporter;

import java.util.List;

import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.Lims;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BaseReportData {

    @NotNull
    public abstract List<PatientTumorLocation> patientTumorLocations();

    @NotNull
    public abstract Lims limsModel();

    @NotNull
    public abstract CenterModel centerModel();

    @NotNull
    public abstract String signaturePath();

   @NotNull
    public abstract String logoPath();
}
