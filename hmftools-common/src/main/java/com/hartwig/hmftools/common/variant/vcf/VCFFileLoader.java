package com.hartwig.hmftools.common.variant.vcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
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
    public static VCFSomaticFile loadSomaticVCF(@NotNull final String basePath, @NotNull final String fileExtension)
            throws IOException, HartwigException {
        return toVCFSomaticFile(loadAllLinesFromVCF(basePath, fileExtension));
    }

    @NotNull
    public static VCFSomaticFile loadSomaticVCF(@NotNull final String file) throws IOException, HartwigException {
        return toVCFSomaticFile(loadAllLinesFromVCF(file));
    }

    @NotNull
    public static List<GermlineVariant> loadGermlineVCF(@NotNull final String basePath,
            @NotNull final String fileExtension) throws IOException, HartwigException {
        final List<String> lines = loadVariantLinesFromVCF(basePath, fileExtension);
        return lines.stream().map(GermlineVariantFactory::fromVCFLine).collect(Collectors.toList());
    }

    @NotNull
    private static VCFSomaticFile toVCFSomaticFile(@NotNull final List<String> lines) {
        final List<String> metaInformationLines = lines.stream().filter(new VCFMetaInformationLinePredicate()).collect(
                Collectors.toList());
        final Optional<String> optHeaderLine = lines.stream().filter(new VCFHeaderLinePredicate()).findFirst();
        assert optHeaderLine.isPresent();
        final String sample = SomaticVariantFactory.sampleFromHeaderLine(optHeaderLine.get());

        final List<String> dataLines = lines.stream().filter(new VCFDataLinePredicate()).collect(Collectors.toList());
        final List<SomaticVariant> variants = dataLines.stream().map(SomaticVariantFactory::fromVCFLine).collect(
                Collectors.toList());
        return new VCFSomaticFile(sample, metaInformationLines, optHeaderLine.get(), variants);
    }

    @NotNull
    private static List<String> loadAllLinesFromVCF(@NotNull final String basePath,
            @NotNull final String fileExtension) throws IOException, HartwigException {
        return FileReader.build().readLines(PathExtensionFinder.build().findPath(basePath, fileExtension));
    }

    @NotNull
    private static List<String> loadAllLinesFromVCF(@NotNull final String file) throws IOException, HartwigException {
        return FileReader.build().readLines(new File(file).toPath());
    }

    @NotNull
    private static List<String> loadVariantLinesFromVCF(@NotNull final String basePath,
            @NotNull final String fileExtension) throws IOException, HartwigException {
        final Path vcfPath = PathExtensionFinder.build().findPath(basePath, fileExtension);
        return LineReader.build().readLines(vcfPath, new VCFDataLinePredicate());
    }
}
