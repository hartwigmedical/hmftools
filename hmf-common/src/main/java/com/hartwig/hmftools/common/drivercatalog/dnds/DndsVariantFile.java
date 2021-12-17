package com.hartwig.hmftools.common.drivercatalog.dnds;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class DndsVariantFile
{
    private static final String DELIMITER = "\t";

    private DndsVariantFile()
    {
    }

    @NotNull
    public static List<DndsVariant> read(@NotNull final String filename) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filename).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull final List<DndsVariant> variants) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(true, variants));
    }

    public static void writeHeader(@NotNull final String filename) throws IOException
    {
        Files.write(new File(filename).toPath(), Collections.singletonList(header()));
    }

    public static void append(@NotNull final String filename, @NotNull final List<DndsVariant> variants) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(false, variants), StandardOpenOption.APPEND);
    }

    @NotNull
    private static List<DndsVariant> fromLines(@NotNull final List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("sample")).map(DndsVariantFile::fromString).collect(Collectors.toList());
    }

    @NotNull
    private static List<String> toLines(boolean header, @NotNull final List<DndsVariant> variants)
    {
        final List<String> lines = Lists.newArrayList();
        if(header)
        {
            lines.add(header());
        }
        variants.stream().map(DndsVariantFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER).add("sample")
                .add("chromosome")
                .add("position")
                .add("ref")
                .add("alt")
                .add("gene")
                .add("worstCodingEffect")
                .add("canonicalCodingEffect")
                .add("biallelic")
                .add("hotspot")
                .add("repeatCount")
                .toString();
    }

    @VisibleForTesting
    @NotNull
    static String toString(@NotNull final DndsVariant variant)
    {
        return new StringJoiner(DELIMITER).add(variant.sampleId())
                .add(variant.chromosome())
                .add(String.valueOf(variant.position()))
                .add(variant.ref())
                .add(variant.alt())
                .add(variant.gene())
                .add(toString(variant.worstCodingEffect()))
                .add(toString(variant.canonicalCodingEffect()))
                .add(String.valueOf(variant.biallelic()))
                .add(String.valueOf(variant.hotspot()))
                .add(String.valueOf(variant.repeatCount()))
                .toString();
    }

    @NotNull
    private static String toString(@NotNull CodingEffect codingEffect)
    {
        return codingEffect != CodingEffect.UNDEFINED ? codingEffect.toString() : Strings.EMPTY;
    }

    @VisibleForTesting
    @NotNull
    static DndsVariant fromString(@NotNull final String line)
    {
        String[] values = line.split(DELIMITER);
        return ImmutableDndsVariant.builder()
                .sampleId(values[0])
                .chromosome(values[1])
                .position(Integer.parseInt(values[2]))
                .ref(values[3])
                .alt(values[4])
                .gene(values[5])
                .worstCodingEffect(toCodingEffect(values[6]))
                .canonicalCodingEffect(toCodingEffect(values[7]))
                .biallelic(Boolean.parseBoolean(values[8]))
                .hotspot(Boolean.parseBoolean(values[9]))
                .repeatCount(Integer.parseInt(values[10]))
                .build();
    }

    @NotNull
    private static CodingEffect toCodingEffect(@NotNull String codingEffect)
    {
        return codingEffect.isEmpty() ? CodingEffect.UNDEFINED : CodingEffect.valueOf(codingEffect);
    }
}
