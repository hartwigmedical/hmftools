package com.hartwig.hmftools.patientreporter;

import java.util.List;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

public final class PatientReport {

    @NotNull
    private final String sample;
    @NotNull
    private final List<VariantReport> variants;
    @NotNull
    private final List<CopyNumberReport> copyNumbers;
    private final int mutationalLoad;

    public PatientReport(@NotNull final String sample, @NotNull final List<VariantReport> variants,
            @NotNull final List<CopyNumberReport> copyNumbers, final int mutationalLoad) {
        this.sample = sample;
        this.variants = variants;
        this.copyNumbers = copyNumbers;
        this.mutationalLoad = mutationalLoad;
    }

    @NotNull
    public String sample() {
        return sample;
    }

    @NotNull
    public List<VariantReport> variants() {
        return variants;
    }

    @NotNull
    public List<CopyNumberReport> copyNumbers() {
        return copyNumbers;
    }

    public int mutationalLoad() {
        return mutationalLoad;
    }
}
