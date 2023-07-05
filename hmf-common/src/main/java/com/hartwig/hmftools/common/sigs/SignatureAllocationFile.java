package com.hartwig.hmftools.common.sigs;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class SignatureAllocationFile
{
    private static final String FILE_EXTENSION = ".sig.allocation.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + FILE_EXTENSION;
    }

    public static List<SignatureAllocation> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<SignatureAllocation> sigAllocations) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(sigAllocations, true));
    }

    public static List<String> toLines(final List<SignatureAllocation> sigAllocations, boolean writerHeader)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        sigAllocations.stream().map(SignatureAllocationFile::toString).forEach(lines::add);
        return lines;
    }

    private static final String SIGNATURE_FLD = "signature";
    private static final String ALLOCATION_FLD = "allocation";
    private static final String PERCENT_FLD = "percent";

    public static List<SignatureAllocation> fromLines(final List<String> lines)
    {
        List<SignatureAllocation> allocations = Lists.newArrayList();

        final String header = lines.get(0);
        lines.remove(0);

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            allocations.add(ImmutableSignatureAllocation.builder()
                    .signature(values[fieldsIndexMap.get(SIGNATURE_FLD)])
                    .allocation(Double.parseDouble(values[fieldsIndexMap.get(ALLOCATION_FLD)]))
                    .percent(Double.parseDouble(values[fieldsIndexMap.get(PERCENT_FLD)]))
                    .build());
        }

        return allocations;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(SIGNATURE_FLD)
                .add(ALLOCATION_FLD)
                .add(PERCENT_FLD)
                .toString();
    }

    private static String toString(final SignatureAllocation sigAllocation)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(sigAllocation.signature()))
                .add(String.format("%.3f", sigAllocation.allocation()))
                .add(String.format("%.5f", sigAllocation.percent()))
                .toString();
    }
}
