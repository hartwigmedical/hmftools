package com.hartwig.hmftools.patientdb.clinical.lims.hospital;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class HospitalPersons {

    @NotNull
    public abstract String hospitalPI();

    // For studies, no requester name is provided.
    @Nullable
    public abstract String requesterName();

    // For studies, no requester email is provided.
    @Nullable
    public abstract String requesterEmail();

}
