package com.hartwig.hmftools.common.purple.purity;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

import org.jetbrains.annotations.NotNull;

public final class PurityContextFile {

    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");
    private static final String DELIMITER = "\t";

    private static final String EXTENSION = ".purple.purity.tsv";
    private static final String EXTENSION_OLD = ".purple.purity";

    @NotNull
    public static PurityQCContext read(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        return readWithQC(PurpleQCFile.generateFilename(basePath, sample), generateFilenameForReading(basePath, sample));
    }

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample) {
        String filename = basePath + File.separator + sample + EXTENSION;
        return (new File(filename).exists()) ? filename : basePath + File.separator + sample + EXTENSION_OLD;
    }

    @NotNull
    public static PurityQCContext readWithQC(@NotNull final String qcFilePath, @NotNull String filePath) throws IOException {
        PurpleQC qc = PurpleQCFile.read(qcFilePath);
        PurityContext purityContext = fromLine(Files.readAllLines(new File(filePath).toPath()).get(1)).build();
        return ImmutablePurityQCContext.builder().purityContext(purityContext).qc(qc).build();
    }

    @NotNull
    public static PurityContext readWithoutQC(@NotNull String filePath) throws IOException {
        return fromLine(Files.readAllLines(new File(filePath).toPath()).get(1)).build();
    }

    @NotNull
    @VisibleForTesting
    static PurityQCContext fromLines(@NotNull List<String> qcLines, @NotNull List<String> fitLines)  {
        final PurpleQC qc = PurpleQCFile.fromLines(qcLines);
        PurityContext purityContext = fromLine(fitLines.get(1)).build();
        return ImmutablePurityQCContext.builder().purityContext(purityContext).qc(qc).build();
    }

    @NotNull
    private static ImmutablePurityContext.Builder fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        ImmutablePurityContext.Builder builder = ImmutablePurityContext.builder()
                .score(score(values))
                .bestFit(bestFit(values))
                .gender(gender(values))
                .method(method(values))
                .polyClonalProportion(polyClonalProportion(values))
                .version(values[14])
                .wholeGenomeDuplication(Boolean.parseBoolean(values[16]))
                .microsatelliteIndelsPerMb(Double.parseDouble(values[17]))
                .microsatelliteStatus(MicrosatelliteStatus.valueOf(values[18]))
                .tumorMutationalLoad(tumorMutationalLoad(values[19]))
                .tumorMutationalLoadStatus(TumorMutationalStatus.valueOf(values[20]))
                .tumorMutationalBurdenPerMb(Double.parseDouble(values[21]))
                .tumorMutationalBurdenStatus(TumorMutationalStatus.valueOf(values[22]))
                .svTumorMutationalBurden(0);

        if (values.length == 24) {
            builder.svTumorMutationalBurden(Integer.parseInt(values[23]));
        }

        return builder;
    }

    public static void write(@NotNull final String basePath, @NotNull final String sample, @NotNull final PurityQCContext context)
            throws IOException {
        PurpleQCFile.write(PurpleQCFile.generateFilename(basePath, sample), context.qc());
        writeBestPurity(basePath, sample, context.purityContext());
    }

    private static void writeBestPurity(@NotNull final String basePath, @NotNull final String sample, @NotNull final PurityContext context)
            throws IOException {
        final String filePath = generateFilenameForWriting(basePath, sample);
        Files.write(new File(filePath).toPath(), toLines(context));
    }

    @NotNull
    private static String generateFilenameForWriting(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull final PurityContext context) {
        return Lists.newArrayList(header(), toString(context));
    }

    @NotNull
    static String header() {
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
                .add("tml")
                .add("tmlStatus")
                .add("tmbPerMb")
                .add("tmbStatus")
                .add("svTumorMutationalBurden")
                .toString();
    }

    @NotNull
    static String toString(@NotNull final PurityContext context) {
        final FittedPurity purity = context.bestFit();
        final FittedPurityScore score = context.score();
        return new StringJoiner(DELIMITER).add(FORMAT.format(purity.purity()))
                .add(FORMAT.format(purity.normFactor()))
                .add(FORMAT.format(purity.score()))
                .add(FORMAT.format(purity.diploidProportion()))
                .add(FORMAT.format(purity.ploidy()))
                .add(String.valueOf(context.gender()))
                .add(String.valueOf(context.method()))
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
                .add(String.valueOf(context.tumorMutationalLoad()))
                .add(String.valueOf(context.tumorMutationalLoadStatus()))
                .add(String.valueOf(context.tumorMutationalBurdenPerMb()))
                .add(String.valueOf(context.tumorMutationalBurdenStatus()))
                .add(String.valueOf(context.svTumorMutationalBurden()))
                .toString();
    }

    @NotNull
    private static FittedPurity bestFit(@NotNull final String[] values) {
        final ImmutableFittedPurity.Builder builder = ImmutableFittedPurity.builder()
                .purity(Double.parseDouble(values[0]))
                .normFactor(Double.parseDouble(values[1]))
                .score(Double.parseDouble(values[2]))
                .diploidProportion(Double.parseDouble(values[3]))
                .ploidy(Double.parseDouble(values[4]))
                .somaticPenalty(0);

        if (values.length > 15) {
            builder.somaticPenalty(Double.parseDouble(values[15]));
        }

        return builder.build();
    }

    @NotNull
    private static Gender gender(@NotNull final String[] values) {
        return Gender.valueOf(values[5]);
    }

    @NotNull
    private static FittedPurityMethod method(@NotNull final String[] values) {
        return FittedPurityMethod.valueOf(values[6]);
    }

    @NotNull
    private static FittedPurityScore score(@NotNull final String[] values) {
        return ImmutableFittedPurityScore.builder()
                .minPurity(Double.parseDouble(values[8]))
                .maxPurity(Double.parseDouble(values[9]))
                .minPloidy(Double.parseDouble(values[10]))
                .maxPloidy(Double.parseDouble(values[11]))
                .minDiploidProportion(Double.parseDouble(values[12]))
                .maxDiploidProportion(Double.parseDouble(values[13]))
                .build();
    }

    private static double polyClonalProportion(@NotNull final String[] values) {
        return Double.parseDouble(values[7]);
    }

    private static int tumorMutationalLoad(@NotNull final String tml) {
        if (tml.contains(".")) {
            return Integer.parseInt(tml.substring(0, tml.indexOf(".")));
        }

        return Integer.parseInt(tml);
    }

}
