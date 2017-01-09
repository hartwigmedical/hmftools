package com.hartwig.hmftools.common.variant.vcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.jetbrains.annotations.NotNull;

public final class VCFFileWriter {

    private static final String VCF_START_LINE = "##fileformat=VCFv4.1";
    // KODU: Variant interpreter figures out the version of the ref genome from below line.
    private static final String VCF_VERSION_LINE = "##reference=GRCh37";
    private static final String VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    public static void writeSomaticVCF(@NotNull final String filePath, @NotNull final List<SomaticVariant> variants)
            throws IOException {
        final Collection<String> vcfLines = Lists.newArrayList();
        vcfLines.addAll(variants.stream().map(SomaticVariantFactory::toVCFLine).collect(Collectors.toList()));
        writeVCF(filePath, vcfLines);
    }

    public static void writeSomaticVCFWithoutInfoFields(@NotNull final String filePath,
            @NotNull final List<SomaticVariant> variants) throws IOException {
        final Collection<String> vcfLines = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            vcfLines.add(SomaticVariantFactory.toVCFLine(variant, true));
        }
        writeVCF(filePath, vcfLines);
    }

    private static void writeVCF(@NotNull final String filePath, @NotNull Collection<String> vcfLines)
            throws IOException {
        final List<String> lines = Lists.newArrayList();
        lines.add(VCF_START_LINE);
        lines.add(VCF_VERSION_LINE);
        lines.add(VCF_HEADER);
        lines.addAll(vcfLines);

        Files.write(new File(filePath).toPath(), lines);
    }
}
