package com.hartwig.hmftools.common.variant.vcfloader;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.io.reader.LineReader;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.jetbrains.annotations.NotNull;

public final class VCFFileLoader {

    private VCFFileLoader() {
    }

    @NotNull
    public static List<SomaticVariant> loadSomaticVCF(@NotNull String basePath, @NotNull String fileExtension)
            throws IOException, HartwigException {
        final List<String> lines = loadVCF(basePath, fileExtension);
        return lines.stream().map(SomaticVariantFactory::fromVCFLine).collect(Collectors.toList());
    }

    @NotNull
    public static List<GermlineVariant> loadGermlineVCF(@NotNull String basePath, @NotNull String fileExtension)
            throws IOException, HartwigException {
        final List<String> lines = loadVCF(basePath, fileExtension);
        return lines.stream().map(GermlineVariantFactory::fromVCFLine).collect(Collectors.toList());
    }

    @NotNull
    private static List<String> loadVCF(@NotNull String basePath, @NotNull String fileExtension)
            throws IOException, HartwigException {
        final Path vcfPath = PathExtensionFinder.build().findPath(basePath, fileExtension);
        return LineReader.build().readLines(vcfPath, new VCFDataLinePredicate());
    }
}
