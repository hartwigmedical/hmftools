package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.cobalt.CobaltCommon.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class MedianRatioFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String EXTENSION = ".cobalt.ratio.median.tsv";

    @NotNull
    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String MEDIAN_RATIO = "medianRatio";
    private static final String COUNT = "count";

    @NotNull
    public static List<MedianRatio> read(final String filename) throws IOException
    {
        List<MedianRatio> ratios = Lists.newArrayList();

        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            ratios.add(ImmutableMedianRatio.builder()
                    .chromosome(values[fieldsIndexMap.get(CHROMOSOME)])
                    .medianRatio(Double.parseDouble(values[fieldsIndexMap.get(MEDIAN_RATIO)]))
                    .count(Integer.parseInt(values[fieldsIndexMap.get(COUNT)])).build());
        }

        return ratios;
    }

    public static void write(final String fileName, List<MedianRatio> ratios) throws IOException
    {
        Files.write(new File(fileName).toPath(), toLines(ratios));
    }

    private static List<String> toLines(final List<MedianRatio> ratio)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        ratio.stream().map(MedianRatioFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "")
                .add(CHROMOSOME)
                .add(MEDIAN_RATIO)
                .add(COUNT).toString();
    }

    private static String toString(final MedianRatio position)
    {
        return new StringJoiner(DELIMITER).add(position.chromosome())
                .add(FORMAT.format(position.medianRatio()))
                .add(String.valueOf(position.count()))
                .toString();
    }
}
