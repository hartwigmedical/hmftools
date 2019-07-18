package com.hartwig.hmftools.common.purple.purity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public final class FittedPurityFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String EXTENSION = ".purple.purity.tsv";
    private static final String EXTENSION_OLD = ".purple.purity";

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + EXTENSION_OLD;
    }

    @NotNull
    public static PurityContext read(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        return read(generateFilenameForReading(basePath, sample));
    }

    @NotNull
    public static PurityContext read(@NotNull String filePath) throws IOException {
        return fromLine(Files.readAllLines(new File(filePath).toPath()).get(1));
    }

    @NotNull
    static PurityContext fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        ImmutablePurityContext.Builder builder = ImmutablePurityContext.builder()
                .score(score(values))
                .bestFit(bestFit(values))
                .gender(gender(values))
                .status(status(values))
                .polyClonalProportion(polyClonalProportion(values))
                .version(values[14])
                .microsatelliteIndelsPerMb(0)
                .microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                .wholeGenomeDuplication(false);

        if (values.length > 16) {
            builder.wholeGenomeDuplication(Boolean.valueOf(values[16]));
        }

        if (values.length > 18) {
            builder.microsatelliteIndelsPerMb(Double.valueOf(values[17]));
            builder.microsatelliteStatus(MicrosatelliteStatus.valueOf(values[18]));
        }

        return builder.build();
    }

    public static void write(@NotNull final String basePath, @NotNull final String sample, @NotNull final PurityContext context)
            throws IOException {
        writeBestPurity(basePath, sample, context);
    }

    private static void writeBestPurity(@NotNull final String basePath, @NotNull final String sample, @NotNull final PurityContext context)
            throws IOException {
        final String filePath = generateFilenameForWriting(basePath, sample);
        Files.write(new File(filePath).toPath(), toLines(context));
    }

    @NotNull
    private static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        //TODO: Once support for reading new / old filename has trickled down to patient report, update this to use new extension!
        return basePath + File.separator + sample + EXTENSION_OLD;
    }

    @NotNull
    static List<String> toLines(@NotNull final PurityContext context) {
        return Lists.newArrayList(header(), toString(context));
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("purity")
                .add("normFactor")
                .add("score")
                .add("diploidProportion")
                .add("ploidy")
                .add("gender")
                .add("status")
                .add("polyclonalProportion")
                .add("minPurity")
                .add("maxPurity")
                .add("minPloidy")
                .add("maxPloidy")
                .add("minDiploidProportion")
                .add("maxDiploidProportion")
                .add("version")
                .add("somaticPenalty")
                .add("wholeGenomeDuplication")
                .add("msIndelsPerMb")
                .add("msStatus")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final PurityContext context) {
        final FittedPurity purity = context.bestFit();
        final FittedPurityScore score = context.score();
        return new StringJoiner(DELIMITER).add(FORMAT.format(purity.purity()))
                .add(FORMAT.format(purity.normFactor()))
                .add(FORMAT.format(purity.score()))
                .add(FORMAT.format(purity.diploidProportion()))
                .add(FORMAT.format(purity.ploidy()))
                .add(String.valueOf(context.gender()))
                .add(String.valueOf(context.status()))
                .add(FORMAT.format(context.polyClonalProportion()))
                .add(FORMAT.format(score.minPurity()))
                .add(FORMAT.format(score.maxPurity()))
                .add(FORMAT.format(score.minPloidy()))
                .add(FORMAT.format(score.maxPloidy()))
                .add(FORMAT.format(score.minDiploidProportion()))
                .add(FORMAT.format(score.maxDiploidProportion()))
                .add(String.valueOf(context.version()))
                .add(FORMAT.format(purity.somaticPenalty()))
                .add(String.valueOf(context.wholeGenomeDuplication()))
                .add(String.valueOf(context.microsatelliteIndelsPerMb()))
                .add(String.valueOf(context.microsatelliteStatus()))
                .toString();
    }

    @NotNull
    private static FittedPurity bestFit(@NotNull final String[] values) {
        final ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder()
                .purity(Double.valueOf(values[0]))
                .normFactor(Double.valueOf(values[1]))
                .score(Double.valueOf(values[2]))
                .diploidProportion(Double.valueOf(values[3]))
                .ploidy(Double.valueOf(values[4]))
                .somaticPenalty(0);

        if (values.length > 15) {
            builder.somaticPenalty(Double.valueOf(values[15]));
        }

        return builder.build();
    }

    @NotNull
    private static Gender gender(@NotNull final String[] values) {
        return Gender.valueOf(values[5]);
    }

    @NotNull
    private static FittedPurityStatus status(@NotNull final String[] values) {
        return FittedPurityStatus.valueOf(values[6]);
    }

    @NotNull
    private static FittedPurityScore score(@NotNull final String[] values) {
        return ImmutableFittedPurityScore.builder()
                .minPurity(Double.valueOf(values[8]))
                .maxPurity(Double.valueOf(values[9]))
                .minPloidy(Double.valueOf(values[10]))
                .maxPloidy(Double.valueOf(values[11]))
                .minDiploidProportion(Double.valueOf(values[12]))
                .maxDiploidProportion(Double.valueOf(values[13]))
                .build();
    }

    private static double polyClonalProportion(@NotNull final String[] values) {
        return Double.valueOf(values[7]);
    }
}
