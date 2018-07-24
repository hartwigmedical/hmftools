package com.hartwig.hmftools.common.variant;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public final class VariantContextFromString {

    private static final VCFCodec CODEC = createTestCodec();

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet("test_sample"));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    private VariantContextFromString() {
    }

    @NotNull
    public static VariantContext decode(@NotNull String line) {
        return createTestCodec().decode(line);
    }
}
