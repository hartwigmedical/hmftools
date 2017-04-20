package com.hartwig.hmftools.common.variant.vcf;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class VCFSomaticFile {

    @NotNull
    private final String sample;
    @NotNull
    private final List<String> originalMetaInformationLines;
    @NotNull
    private final String originalHeaderLine;
    @NotNull
    private final List<SomaticVariant> variants;

    VCFSomaticFile(@NotNull final String sample, @NotNull final List<String> originalMetaInformationLines,
            @NotNull final String originalHeaderLine, @NotNull final List<SomaticVariant> variants) {
        this.sample = sample;
        this.originalMetaInformationLines = originalMetaInformationLines;
        this.originalHeaderLine = originalHeaderLine;
        this.variants = variants;
    }

    @NotNull
    public String sample() {
        return sample;
    }

    @NotNull
    List<String> originalMetaInformationLines() {
        return originalMetaInformationLines;
    }

    @NotNull
    String originalHeaderLine() {
        return originalHeaderLine;
    }

    @NotNull
    public List<SomaticVariant> variants() {
        return variants;
    }
}
