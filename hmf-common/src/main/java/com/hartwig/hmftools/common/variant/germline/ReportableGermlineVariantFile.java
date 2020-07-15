package com.hartwig.hmftools.common.variant.germline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

// TODO as part of patient reporter upgrade to 7.14 (DEV-1055):
//  - Remove 'program' column
//  - Remove 'filter' column -> assume always PASS
//  - Rename 'alts' column to 'alt'
//  - Potentially remove other fields that are not used downstream.
public final class ReportableGermlineVariantFile
{
    private static final String DELIMITER = "\t";

    private static final String FILE_EXTENSION = ".reportable_germline_variant.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull List<ReportableGermlineVariant> dataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(dataList));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<ReportableGermlineVariant> dataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        dataList.stream().map(ReportableGermlineVariantFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("program")
                .add("gene")
                .add("chromosome")
                .add("position")
                .add("ref")
                .add("alts")
                .add("filter")
                .add("codingEffect")
                .add("hgvsCoding")
                .add("hgvsProtein")
                .add("alleleReadCount")
                .add("totalReadCount")
                .add("adjustedVaf")
                .add("adjustedCopyNumber")
                .add("biallelic")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final ReportableGermlineVariant data)
    {
        return new StringJoiner(DELIMITER)
                .add("HMF")
                .add(String.valueOf(data.gene()))
                .add(String.valueOf(data.chromosome()))
                .add(String.valueOf(data.position()))
                .add(String.valueOf(data.ref()))
                .add(String.valueOf(data.alt()))
                .add(String.valueOf(data.passFilter()))
                .add(String.valueOf(data.codingEffect()))
                .add(String.valueOf(data.hgvsCoding()))
                .add(String.valueOf(data.hgvsProtein()))
                .add(String.valueOf(data.alleleReadCount()))
                .add(String.valueOf(data.totalReadCount()))
                .add(String.valueOf(data.adjustedVaf()))
                .add(String.valueOf(data.adjustedCopyNumber()))
                .add(String.valueOf(data.biallelic()))
                .toString();
    }
}
