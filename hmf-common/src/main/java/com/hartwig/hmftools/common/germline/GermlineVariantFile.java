package com.hartwig.hmftools.common.germline;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.jetbrains.annotations.NotNull;

public class GermlineVariantFile
{
    private static final String DELIMITER = "\t";

    private static final String FILE_EXTENSION = ".bachelor.germline_variant.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<GermlineVariant> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<GermlineVariant> dataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(dataList));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<GermlineVariant> dataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        dataList.stream().map(GermlineVariantFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<GermlineVariant> fromLines(@NotNull List<String> lines)
    {
        return lines.stream()
                .skip(1)
                .map(GermlineVariantFile::fromString).collect(toList());
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER)
                .add("chromosome")
                .add("position")
                .add("filter")
                .add("reported")
                .add("pathogenic")
                .add("type")
                .add("ref")
                .add("alts")
                .add("gene")
                .add("transcriptId")
                .add("effects")
                .add("codingEffect")
                .add("microhomology")
                .add("repeatSequence")
                .add("repeatCount")
                .add("alleleReadCount")
                .add("totalReadCount")
                .add("adjustedVaf")
                .add("adjustedCopyNumber")
                .add("trinucleotideContext")
                .add("hgvsProtein")
                .add("hgvsCoding")
                .add("biallelic")
                .add("minorAllelePloidy")
                .add("refStatus")
                .add("codonInfo")
                .add("clinvarInfo")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final GermlineVariant data)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(data.chromosome()))
                .add(String.valueOf(data.position()))
                .add(String.valueOf(data.filter()))
                .add(String.valueOf(data.reported()))
                .add(String.valueOf(data.pathogenic()))
                .add(String.valueOf(data.type()))
                .add(String.valueOf(data.ref()))
                .add(String.valueOf(data.alts()))
                .add(String.valueOf(data.gene()))
                .add(String.valueOf(data.transcriptId()))
                .add(String.valueOf(data.effects()))
                .add(String.valueOf(data.codingEffect()))
                .add(String.valueOf(data.microhomology()))
                .add(String.valueOf(data.repeatSequence()))
                .add(String.valueOf(data.repeatCount()))
                .add(String.valueOf(data.alleleReadCount()))
                .add(String.valueOf(data.totalReadCount()))
                .add(String.valueOf(data.adjustedVaf()))
                .add(String.valueOf(data.adjustedCopyNumber()))
                .add(String.valueOf(data.trinucleotideContext()))
                .add(String.valueOf(data.hgvsProtein()))
                .add(String.valueOf(data.hgvsCoding()))
                .add(String.valueOf(data.biallelic()))
                .add(String.valueOf(data.minorAlleleJcn()))
                .add(String.valueOf(data.refStatus()))
                .add(String.valueOf(data.clinvarInfo()))
                .toString();
    }

    @NotNull
    private static GermlineVariant fromString(@NotNull final String data)
    {
        String[] values = data.split(DELIMITER, -1);

        int index = 0;

        return ImmutableGermlineVariant.builder()
                .chromosome(values[index++])
                .position(Long.parseLong(values[index++]))
                .filter(FilterType.valueOf(values[index++]))
                .reported(Boolean.parseBoolean(values[index++]))
                .pathogenic(PathogenicType.valueOf(values[index++]))
                .type(values[index++])
                .ref(values[index++])
                .alts(values[index++])
                .gene(values[index++])
                .transcriptId(values[index++])
                .effects(values[index++])
                .codingEffect(CodingEffect.valueOf(values[index++]))
                .microhomology(values[index++])
                .repeatSequence(values[index++])
                .repeatCount(Integer.parseInt(values[index++]))
                .alleleReadCount(Integer.parseInt(values[index++]))
                .totalReadCount(Integer.parseInt(values[index++]))
                .adjustedVaf(Double.parseDouble(values[index++]))
                .adjustedCopyNumber(Double.parseDouble(values[index++]))
                .trinucleotideContext(values[index++])
                .hgvsProtein(values[index++])
                .hgvsCoding(values[index++])
                .biallelic(Boolean.parseBoolean(values[index++]))
                .minorAlleleJcn(Double.parseDouble(values[index++]))
                .refStatus(values[index++])
                .clinvarInfo(values[index++])
                .build();
    }
}
