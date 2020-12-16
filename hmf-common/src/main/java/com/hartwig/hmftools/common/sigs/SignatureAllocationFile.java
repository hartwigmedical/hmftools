package com.hartwig.hmftools.common.sigs;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class SignatureAllocationFile
{
    public static final String DELIMITER = "\t";

    private static final String FILE_EXTENSION = ".sig.allocation.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<SignatureAllocation> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<SignatureAllocation> sigAllocations) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(sigAllocations, true));
    }

    @NotNull
    public static List<String> toLines(@NotNull final List<SignatureAllocation> sigAllocations, boolean writerHeader)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        sigAllocations.stream().map(SignatureAllocationFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    public static List<SignatureAllocation> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("signature")).map(SignatureAllocationFile::fromString).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("signature")
                .add("allocation")
                .add("percent")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final SignatureAllocation sigAllocation)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(sigAllocation.signature()))
                .add(String.format("%.3f", sigAllocation.allocation()))
                .add(String.format("%.5f", sigAllocation.percent()))
                .toString();
    }

    @NotNull
    private static SignatureAllocation fromString(@NotNull final String sigAllocation)
    {
        String[] values = sigAllocation.split(DELIMITER);

        int index = 0;

        return ImmutableSignatureAllocation.builder()
                .signature(values[index++])
                .allocation(Double.parseDouble(values[index++]))
                .percent(Double.parseDouble(values[index++]))
                .build();
    }
}
