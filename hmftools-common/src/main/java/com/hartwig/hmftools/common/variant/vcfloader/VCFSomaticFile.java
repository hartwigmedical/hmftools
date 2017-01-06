package com.hartwig.hmftools.common.variant.vcfloader;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class VCFSomaticFile {

    @NotNull
    private final String sample;
    @NotNull
    private final List<SomaticVariant> variants;

    VCFSomaticFile(@NotNull final String sample, @NotNull final List<SomaticVariant> variants) {
        this.sample = sample;
        this.variants = variants;
    }

    @NotNull
    public String sample() {
        return sample;
    }

    @NotNull
    public List<SomaticVariant> variants() {
        return variants;
    }
}
