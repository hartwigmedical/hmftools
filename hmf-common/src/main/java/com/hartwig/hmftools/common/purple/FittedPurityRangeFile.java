package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.purple.PurpleCommon.DELIMITER;
import static com.hartwig.hmftools.common.purple.BestFit.bestFitPerPurity;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class FittedPurityRangeFile
{
    private static final DecimalFormat FORMAT = PurpleCommon.decimalFormat("0.0000");

    private static final String EXTENSION = ".purple.purity.range.tsv";

    @NotNull
    public static String generateFilenameForReading(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    private static final String PURITY = "purity";
    private static final String NORM_FACTOR = "normFactor";
    private static final String SCORE = "score";
    private static final String DIPLOD_PROPORTION = "diploidProportion";
    private static final String PLOIDY = "ploidy";
    private static final String SOMATIC_PENALTY = "somaticPenalty";

    @NotNull
    public static List<FittedPurity> readAll(final String basePath, final String sample) throws IOException
    {
        final String filePath = generateFilenameForReading(basePath, sample);

        List<FittedPurity> fittedPurities = Lists.newArrayList();
        List<String> lines = Files.readAllLines(new File(filePath).toPath());

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            fittedPurities.add(ImmutableFittedPurity.builder()
                    .purity(Double.parseDouble(values[fieldsIndexMap.get(PURITY)]))
                    .normFactor(Double.parseDouble(values[fieldsIndexMap.get(NORM_FACTOR)]))
                    .score(Double.parseDouble(values[fieldsIndexMap.get(SCORE)]))
                    .diploidProportion(Double.parseDouble(values[fieldsIndexMap.get(DIPLOD_PROPORTION)]))
                    .ploidy(Double.parseDouble(values[fieldsIndexMap.get(PLOIDY)]))
                    .somaticPenalty(Double.parseDouble(values[fieldsIndexMap.get(SOMATIC_PENALTY)]))
                    .build());
        }

        return fittedPurities;
    }

    @NotNull
    public static List<FittedPurity> readBestFitPerPurity(final String basePath, final String sample) throws IOException
    {
        return bestFitPerPurity(readAll(basePath, sample));
    }

    public static void write(final String basePath, final String sample, final List<FittedPurity> purity)
            throws IOException
    {
        final String filePath = generateFilenameForWriting(basePath, sample);
        Files.write(new File(filePath).toPath(), toLines(purity));
    }

    @NotNull
    private static String generateFilenameForWriting(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    @VisibleForTesting
    private static List<String> toLines(final List<FittedPurity> purity)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(FittedPurityRangeFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "")
                .add(PURITY)
                .add(NORM_FACTOR)
                .add(SCORE)
                .add(DIPLOD_PROPORTION)
                .add(PLOIDY)
                .add(SOMATIC_PENALTY)
                .toString();
    }

    @NotNull
    private static String toString(final FittedPurity purity)
    {
        return new StringJoiner(DELIMITER)
                .add(FORMAT.format(purity.purity()))
                .add(FORMAT.format(purity.normFactor()))
                .add(FORMAT.format(purity.score()))
                .add(FORMAT.format(purity.diploidProportion()))
                .add(FORMAT.format(purity.ploidy()))
                .add(FORMAT.format(purity.somaticPenalty()))
                .toString();
    }
}
