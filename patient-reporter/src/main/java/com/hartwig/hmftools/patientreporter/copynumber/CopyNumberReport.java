package com.hartwig.hmftools.patientreporter.copynumber;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.gene.GeneRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CopyNumberReport implements GeneRegion, CopyNumber {

    @NotNull
    public abstract CopyNumberReportType type();

    @NotNull
    public String description() {
        return type().description();
    }

    @NotNull
    public String position() {
        return chromosome() + chromosomeBand();
    }

}
