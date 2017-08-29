package com.hartwig.hmftools.common.variant.vcf;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Function;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
import com.hartwig.hmftools.common.io.reader.FileReader;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.GermlineVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticTruthSetVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariant;
import com.hartwig.hmftools.common.variant.strelka.StrelkaSomaticVariantFactory;

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
    public static VCFSomaticTruthSetFile loadSomaticTruthSetVCF(@NotNull final String file) throws IOException, HartwigException {
        return toVCFSomaticTruthSet(loadAllLinesFromVCF(file));
    }

    @NotNull
    public static StrelkaVCFSomaticFile loadStrelkaVCF(@NotNull final String file) throws IOException, HartwigException {
        return toStrelkaVCF(loadAllLinesFromVCF(file));
    }

    @NotNull
    public static VCFGermlineFile loadGermlineVCF(@NotNull final String file) throws IOException, HartwigException {
        return toVCFGermlineFile(loadAllLinesFromVCF(file));
    }

    @NotNull
    private static VCFSomaticFile toVCFSomaticFile(@NotNull final List<String> lines) {
        final List<String> metaInformationLines = extractMetaInformation(lines);
        final String header = extractHeader(lines);
        final String sample = SomaticVariantFactory.sampleFromHeaderLine(header);
        final List<SomaticVariant> variants = variants(lines, SomaticVariantFactory::fromVCFLine);
        return ImmutableVCFSomaticFile.builder()
                .sample(sample)
                .originalMetaInformationLines(metaInformationLines)
                .originalHeaderLine(header)
                .variants(variants)
                .build();
    }

    @NotNull
    private static VCFSomaticTruthSetFile toVCFSomaticTruthSet(@NotNull final List<String> lines) {
        final List<String> metaInformationLines = extractMetaInformation(lines);
        final String header = extractHeader(lines);
        final String sample = SomaticVariantFactory.sampleFromHeaderLine(header);
        final List<SomaticTruthSetVariant> variants = variants(lines, SomaticVariantFactory::fromTruthSetVCFLine);
        return ImmutableVCFSomaticTruthSetFile.builder()
                .sample(sample)
                .originalMetaInformationLines(metaInformationLines)
                .originalHeaderLine(header)
                .variants(variants)
                .build();
    }

    @NotNull
    private static StrelkaVCFSomaticFile toStrelkaVCF(@NotNull final List<String> lines) {
        final List<String> metaInformationLines = extractMetaInformation(lines);
        final String header = extractHeader(lines);
        final List<StrelkaSomaticVariant> variants = variants(lines, StrelkaSomaticVariantFactory::fromVCFLine);
        return ImmutableStrelkaVCFSomaticFile.builder()
                .originalMetaInformationLines(metaInformationLines)
                .originalHeaderLine(header)
                .variants(variants)
                .build();
    }

    @NotNull
    private static VCFGermlineFile toVCFGermlineFile(@NotNull final List<String> lines) {
        final List<String> metaInformationLines = extractMetaInformation(lines);
        final String header = extractHeader(lines);
        final String refSample = GermlineVariantFactory.refSampleFromHeaderLine(header);
        final String tumorSample = GermlineVariantFactory.tumorSampleFromHeaderLine(header);
        final List<GermlineVariant> variants = variants(lines, GermlineVariantFactory::fromVCFLine);
        return ImmutableVCFGermlineFile.builder()
                .refSample(refSample)
                .tumorSample(tumorSample)
                .originalMetaInformationLines(metaInformationLines)
                .originalHeaderLine(header)
                .variants(variants)
                .build();
    }

    @NotNull
    private static List<String> loadAllLinesFromVCF(@NotNull final String basePath, @NotNull final String fileExtension)
            throws IOException, HartwigException {
        return FileReader.build().readLines(PathExtensionFinder.build().findPath(basePath, fileExtension));
    }

    @NotNull
    private static List<String> loadAllLinesFromVCF(@NotNull final String file) throws IOException, HartwigException {
        return FileReader.build().readLines(new File(file).toPath());
    }

    @NotNull
    private static String extractHeader(@NotNull final List<String> lines) {
        final Optional<String> optHeaderLine = lines.stream().filter(new VCFHeaderLinePredicate()).findFirst();
        Preconditions.checkState(optHeaderLine.isPresent());
        return optHeaderLine.get();
    }

    @NotNull
    private static List<String> extractMetaInformation(@NotNull final List<String> lines) {
        return lines.stream().filter(new VCFMetaInformationLinePredicate()).collect(toList());
    }

    @NotNull
    private static <T extends Variant> List<T> variants(@NotNull final List<String> lines, Function<String, T> transform) {
        return lines.stream().filter(new VCFDataLinePredicate()).map(transform).filter(Objects::nonNull).collect(toList());
    }
}
