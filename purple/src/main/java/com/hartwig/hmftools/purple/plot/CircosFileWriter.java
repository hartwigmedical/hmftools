package com.hartwig.hmftools.purple.plot;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.StringJoiner;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

import org.jetbrains.annotations.Nullable;

public final class CircosFileWriter
{
    public static <T extends GenomeRegion> void writeRegions(
            final String filePath, Collection<T> values, ToDoubleFunction<T> valueExtractor) throws IOException
    {
        Function<T, String> toString = t -> transformRegion(valueExtractor, null, t);
        writeCircosFile(filePath, values, toString);
    }

    public static <T extends GenomeRegion> void writeRegions(
            final String filePath,
            Collection<T> values,
            ToDoubleFunction<T> valueExtractor,
            ToDoubleFunction<T> glyphSizeExtractor) throws IOException
    {
        Function<T, String> toString = t -> transformRegion(valueExtractor, glyphSizeExtractor, t);
        writeCircosFile(filePath, values, toString);
    }

    public static <T extends GenomePosition> void writePositions(
            final String filePath, Collection<T> values, ToDoubleFunction<T> valueExtractor) throws IOException
    {
        Function<T, String> toString = t -> transformPosition(valueExtractor, t);
        writeCircosFile(filePath, values, toString);
    }

    public static <T> void writeCircosFile(
            final String filePath, Collection<T> values, Function<T, String> toStringFunction) throws IOException
    {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        values.stream().map(toStringFunction).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    private static String header()
    {
        return "#chromosome\tstart\tend\tvalue";
    }

    private static <T extends GenomePosition> String transformPosition(ToDoubleFunction<T> valueExtractor, T position)
    {
        return new StringBuilder().append(circosContig(position.chromosome()))
                .append('\t')
                .append(position.position())
                .append('\t')
                .append(position.position())
                .append('\t')
                .append(valueExtractor.applyAsDouble(position))
                .toString();
    }

    private static <T extends GenomeRegion> String transformRegion(
            ToDoubleFunction<T> valueExtractor,
            @Nullable ToDoubleFunction<T> glyphSizeExtractor, T region)
    {
        StringBuilder stringBuilder = new StringBuilder().append(circosContig(region.chromosome()))
                .append('\t')
                .append(region.start())
                .append('\t')
                .append(region.end())
                .append('\t')
                .append(String.format("%.4f", valueExtractor.applyAsDouble(region)));

        if(glyphSizeExtractor != null)
        {
            stringBuilder.append('\t')
                    .append(String.format("glyph_size=%.2f", glyphSizeExtractor.applyAsDouble(region)));
        }
        return stringBuilder.toString();
    }

    public static String circosContig(String chromosome)
    {
        return "hs" + HumanChromosome.fromString(chromosome);
    }

    static String transformPosition(VariantContextDecorator position, Function<VariantContextDecorator, String> colourFunction)
    {
        return new StringJoiner("\t").add(circosContig(position.chromosome()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.adjustedVaf()))
                .add("color=" + colourFunction.apply(position))
                .toString();
    }
}
