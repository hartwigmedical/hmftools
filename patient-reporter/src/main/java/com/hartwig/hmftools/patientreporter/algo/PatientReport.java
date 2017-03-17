package com.hartwig.hmftools.patientreporter.algo;

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
    @NotNull
    private final String tumorType;
    private final double tumorPercentage;

    public PatientReport(@NotNull final String sample, @NotNull final List<VariantReport> variants,
            @NotNull final List<CopyNumberReport> copyNumbers, final int mutationalLoad,
            @NotNull final String tumorType, final double tumorPercentage) {
        this.sample = sample;
        this.variants = variants;
        this.copyNumbers = copyNumbers;
        this.mutationalLoad = mutationalLoad;
        this.tumorType = tumorType;
        this.tumorPercentage = tumorPercentage;
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

    @NotNull
    public String tumorType() {
        return tumorType;
    }

    public String tumorPercentageString() {
        return Long.toString(Math.round(tumorPercentage * 100D)) + "%";
    }
}
