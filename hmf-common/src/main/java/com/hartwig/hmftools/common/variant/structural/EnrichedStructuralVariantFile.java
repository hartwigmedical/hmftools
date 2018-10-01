package com.hartwig.hmftools.common.variant.structural;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

import static java.util.stream.Collectors.toList;

public enum EnrichedStructuralVariantFile {
    ;

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".sv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static List<EnrichedStructuralVariant> read(final String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull Collection<EnrichedStructuralVariant> copyNumbers) throws IOException {
        Files.write(new File(filename).toPath(), toLines(copyNumbers));
    }

    @NotNull
    static List<String> toLines(@NotNull final Collection<EnrichedStructuralVariant> variants) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        variants.stream().map(EnrichedStructuralVariantFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<EnrichedStructuralVariant> fromLines(@NotNull Collection<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(EnrichedStructuralVariantFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "")
                .add("id")
                .add("startChromosome")
                .add("endChromosome")
                .add("startPosition")
                .add("endPosition")
                .add("startOrientation")
                .add("endOrientation")
                .add("startHomologySequence")
                .add("endHomologySequence")
                .add("startAF")
                .add("endAF")
                .add("ploidy")
                .add("adjustedStartAF")
                .add("adjustedEndAF")
                .add("adjustedStartCopyNumber")
                .add("adjustedEndCopyNumber")
                .add("adjustedStartCopyNumberChange")
                .add("adjustedEndCopyNumberChange")
                .add("insertSequence")
                .add("type")
                .add("filter")
                .add("imprecise")
                .add("qualScore")
                .add("event")
                .add("startTumourVariantFragmentCount")
                .add("startTumourReferenceFragmentCount")
                .add("startNormalVariantFragmentCount")
                .add("startNormalReferenceFragmentCount")
                .add("endTumourVariantFragmentCount")
                .add("endTumourReferenceFragmentCount")
                .add("endNormalVariantFragmentCount")
                .add("endNormalReferenceFragmentCount")
                .add("startIntervalOffsetStart")
                .add("startIntervalOffsetEnd")
                .add("endIntervalOffsetStart")
                .add("endIntervalOffsetEnd")
                .add("inexactHomologyOffsetStart")
                .add("inexactHomologyOffsetEnd")
                .add("startLinkedBy")
                .add("endLinkedBy")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final EnrichedStructuralVariant variant) {
        return new StringJoiner(DELIMITER)
                .add(variant.id())
                .add(variant.start().chromosome())
                .add(variant.end() == null ? null : variant.end().chromosome())
                .add(String.valueOf(variant.start().position()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().position()))
                .add(Integer.toString(variant.start().orientation()))
                .add(variant.end() == null ? null : Integer.toString(variant.end().orientation()))
                .add(variant.start().homology())
                .add(variant.end() == null ? null : variant.end().homology())
                .add(String.valueOf(variant.start().alleleFrequency()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().alleleFrequency()))
                .add(String.valueOf(variant.ploidy()))
                .add(String.valueOf(variant.start().adjustedAlleleFrequency()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().adjustedAlleleFrequency()))
                .add(String.valueOf(variant.start().adjustedCopyNumber()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().adjustedCopyNumber()))
                .add(String.valueOf(variant.start().adjustedCopyNumberChange()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().adjustedCopyNumberChange()))
                .add(String.valueOf(variant.insertSequence()))
                .add(String.valueOf(variant.type()))
                .add(String.valueOf(variant.filter()))
                .add(String.valueOf(variant.imprecise()))
                .add(String.valueOf(variant.qualityScore()))
                .add(String.valueOf(variant.event()))
                .add(String.valueOf(variant.start().tumourVariantFragmentCount()))
                .add(String.valueOf(variant.start().tumourReferenceFragmentCount()))
                .add(String.valueOf(variant.start().normalVariantFragmentCount()))
                .add(String.valueOf(variant.start().normalReferenceFragmentCount()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().tumourVariantFragmentCount()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().tumourReferenceFragmentCount()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().normalVariantFragmentCount()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().normalReferenceFragmentCount()))
                .add(String.valueOf(variant.start().startOffset()))
                .add(String.valueOf(variant.start().endOffset()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().startOffset()))
                .add(String.valueOf(variant.end() == null ? null : variant.end().endOffset()))
                .add(String.valueOf(variant.start().inexactHomologyOffsetStart()))
                .add(String.valueOf(variant.start().inexactHomologyOffsetEnd()))
                .add(variant.startLinkedBy())
                .add(variant.endLinkedBy())
                .toString();
    }

    @NotNull
    private static EnrichedStructuralVariant fromString(@NotNull final String variant) {
        throw new RuntimeException("Not Yet Implemented");
    }
}
