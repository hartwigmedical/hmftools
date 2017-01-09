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
import org.jetbrains.annotations.Nullable;

public final class VCFFileWriter {

    private static final String VCF_START_LINE = "##fileformat=VCFv4.1";
    // KODU: Variant interpreter figures out the version of the ref genome from below line.
    private static final String VCF_VERSION_LINE = "##reference=GRCh37";
    private static final String VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    private static final String VCF_INFO_FIELD = "CLEAN";
    private static final String VCF_OPTIONAL_HEADER = "##INFO=<ID=" + VCF_INFO_FIELD
            + ",Number=1,Type=Integer,Description=\"Indicates original INFO field has been cleaned\">";

    private static final String VCF_FORMAT_GT = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
    private static final String VCF_FORMAT_AD = "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">";
    private static final String VCF_FORMAT_DP = "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">";

    public static void writeSomaticVCF(@NotNull final String filePath, @NotNull final List<SomaticVariant> variants)
            throws IOException {
        final Collection<String> vcfLines = Lists.newArrayList();
        vcfLines.addAll(variants.stream().map(SomaticVariantFactory::toVCFLine).collect(Collectors.toList()));
        writeVCF(filePath, vcfLines, null);
    }

    // KODU: This function is supposed to generate a VCF which can be merged with other similar VCFs using GATK
    // CombineVariants, but so far it does not really work well!
    public static void writeCleanedSomaticVCF(@NotNull final String filePath,
            @NotNull final List<SomaticVariant> variants) throws IOException {
        final Collection<String> vcfLines = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            vcfLines.add(SomaticVariantFactory.toVCFLine(variant, VCF_INFO_FIELD + "=1"));
        }
        writeVCF(filePath, vcfLines, VCF_OPTIONAL_HEADER);
    }

    private static void writeVCF(@NotNull final String filePath, @NotNull Collection<String> vcfLines,
            @Nullable String optionalHeaderLine) throws IOException {
        final List<String> lines = Lists.newArrayList();
        lines.add(VCF_START_LINE);
        lines.add(VCF_VERSION_LINE);

        if (optionalHeaderLine != null) {
            lines.add(optionalHeaderLine);
        }

        lines.add(VCF_FORMAT_GT);
        lines.add(VCF_FORMAT_AD);
        lines.add(VCF_FORMAT_DP);

        lines.add(VCF_HEADER);

        lines.addAll(vcfLines);

        Files.write(new File(filePath).toPath(), lines);
    }
}
