package com.hartwig.hmftools.common.test;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public final class VariantContextFromString {

    public static final String SAMPLE = "test_sample";

    private VariantContextFromString() {
    }

    @NotNull
    public static VariantContext decode(@NotNull String line) {
        return decode(SAMPLE, line);
    }

    @NotNull
    public static VariantContext decode(@NotNull String sample, @NotNull String line) {
        return createTestCodec(sample).decode(line);
    }

    @NotNull
    private static VCFCodec createTestCodec(@NotNull String sample) {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(sample));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
