package com.hartwig.hmftools.common.purple.gene;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String EXTENSION = ".purple.cnv.gene.tsv";

    private GeneCopyNumberFile() {
    }

    @NotNull
    public static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    public static List<GeneCopyNumber> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String fileName, @NotNull List<GeneCopyNumber> geneCopyNumbers) throws IOException {
        Files.write(new File(fileName).toPath(), toLines(geneCopyNumbers));
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final List<GeneCopyNumber> ratio) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(GeneCopyNumberFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "")
                .add("chromosome")
                .add("start")
                .add("end")
                .add("gene")
                .add("minCopyNumber")
                .add("maxCopyNumber")
                .add("unused")
                .add("somaticRegions")
                .add("germlineHomDeletionRegions")
                .add("germlineHetToHomDeletionRegions")
                .add("transcriptId")
                .add("chromosomeBand")
                .add("minRegions")
                .add("minRegionStart")
                .add("minRegionEnd")
                .add("minRegionStartSupport")
                .add("minRegionEndSupport")
                .add("minRegionMethod")
                .add("minMinorAlleleCopyNumber")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final GeneCopyNumber geneCopyNumber) {

        return new StringJoiner(DELIMITER).add(geneCopyNumber.chromosome())
                .add(String.valueOf(geneCopyNumber.start()))
                .add(String.valueOf(geneCopyNumber.end()))
                .add(geneCopyNumber.gene())
                .add(FORMAT.format(geneCopyNumber.minCopyNumber()))
                .add(FORMAT.format(geneCopyNumber.maxCopyNumber()))
                .add(String.valueOf(0))
                .add(String.valueOf(geneCopyNumber.somaticRegions()))
                .add(String.valueOf(geneCopyNumber.germlineHomRegions()))
                .add(String.valueOf(geneCopyNumber.germlineHet2HomRegions()))
                .add(geneCopyNumber.transcriptID())
                .add(geneCopyNumber.chromosomeBand())
                .add(String.valueOf(geneCopyNumber.minRegions()))
                .add(String.valueOf(geneCopyNumber.minRegionStart()))
                .add(String.valueOf(geneCopyNumber.minRegionEnd()))
                .add(String.valueOf(geneCopyNumber.minRegionStartSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionEndSupport()))
                .add(String.valueOf(geneCopyNumber.minRegionMethod()))
                .add(FORMAT.format(geneCopyNumber.minMinorAlleleCopyNumber()))
                .toString();
    }

    @NotNull
    @VisibleForTesting
    static List<GeneCopyNumber> fromLines(@NotNull List<String> lines) {

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
        lines.remove(0);

        int chrIndex = fieldsIndexMap.get("chromosome");
        int startIndex = fieldsIndexMap.get("start");
        int endIndex = fieldsIndexMap.get("end");
        int geneIndex = fieldsIndexMap.get("gene");
        int minCnIndex = fieldsIndexMap.get("minCopyNumber");
        int maxCnIndex = fieldsIndexMap.get("maxCopyNumber");
        int somRegionsIndex = fieldsIndexMap.get("somaticRegions");
        int germHmDelRegionsIndex = fieldsIndexMap.get("germlineHomDeletionRegions");
        int germHtHDelRegionsIndex = fieldsIndexMap.get("germlineHetToHomDeletionRegions");
        int transIdIndex = fieldsIndexMap.get("transcriptId");
        int chrBandIndex = fieldsIndexMap.get("chromosomeBand");
        int minRegionIndex = fieldsIndexMap.get("minRegions");
        int minRegionStartIndex = fieldsIndexMap.get("minRegionStart");
        int minRegionEndIndex = fieldsIndexMap.get("minRegionEnd");
        int minRegionStartSupIndex = fieldsIndexMap.get("minRegionStartSupport");
        int minRegionEndSupIndex = fieldsIndexMap.get("minRegionEndSupport");
        int minRegionMethodIndex = fieldsIndexMap.get("minRegionMethod");
        int mmACnIndex = fieldsIndexMap.get("minMinorAlleleCopyNumber");

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();

        for(final String line : lines) {
            String[] values = line.split(DELIMITER, -1);

            final ImmutableGeneCopyNumber.Builder builder = ImmutableGeneCopyNumber.builder()
                    .chromosome(values[chrIndex])
                    .start(Long.parseLong(values[startIndex]))
                    .end(Long.parseLong(values[endIndex]))
                    .gene(values[geneIndex])
                    .minCopyNumber(Double.parseDouble(values[minCnIndex]))
                    .maxCopyNumber(Double.parseDouble(values[maxCnIndex]))
                    .somaticRegions(Integer.parseInt(values[somRegionsIndex]))
                    .germlineHomRegions(Integer.parseInt(values[germHmDelRegionsIndex]))
                    .germlineHet2HomRegions(Integer.parseInt(values[germHtHDelRegionsIndex]))
                    .transcriptID(values[transIdIndex])
                    .chromosomeBand(values[chrBandIndex])
                    .minRegions(minRegionIndex)
                    .minRegionStart(Long.parseLong(values[minRegionStartIndex]))
                    .minRegionEnd(Long.parseLong(values[minRegionEndIndex]))
                    .minRegionMethod(CopyNumberMethod.valueOf(values[minRegionMethodIndex]))
                    .minRegionStartSupport(SegmentSupport.valueOf(values[minRegionStartSupIndex]))
                    .minRegionEndSupport(SegmentSupport.valueOf(values[minRegionEndSupIndex]))
                    .minMinorAlleleCopyNumber(Double.parseDouble(values[mmACnIndex]));

            geneCopyNumbers.add(builder.build());
        }

        return geneCopyNumbers;
    }
}
