package com.hartwig.hmftools.common.teal;

import java.io.File;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.commons.lang3.Validate;

public class TelomereLengthFile
{
    enum Column
    {
        sampleId, type, rawTelomereLength, finalTelomereLength, fullFragments, cRichPartialFragments, gRichPartialFragments,
        totalTelomericReads, purity, ploidy, duplicateProportion, meanReadsPerKb, gc50ReadsPerKb
    }

    private static final String FILE_EXTENSION = ".teal.tellength.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static TelomereLength read(final String filename)
    {
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            // there should only be one row
            List<TelomereLength> telomereLengths = reader.stream().map(row -> ImmutableTelomereLength.builder()
                    .type(Objects.requireNonNull(row.get(Column.type)))
                    .rawTelomereLength(Objects.requireNonNull(row.getDouble(Column.rawTelomereLength)))
                    .finalTelomereLength(Objects.requireNonNull(row.getDouble(Column.finalTelomereLength)))
                    .fullFragments(Objects.requireNonNull(row.getInt(Column.fullFragments)))
                    .cRichPartialFragments(Objects.requireNonNull(row.getInt(Column.cRichPartialFragments)))
                    .gRichPartialFragments(Objects.requireNonNull(row.getInt(Column.gRichPartialFragments)))
                    .totalTelomericReads(Objects.requireNonNull(row.getInt(Column.totalTelomericReads)))
                    .purity(Objects.requireNonNull(row.getDouble(Column.purity)))
                    .ploidy(Objects.requireNonNull(row.getDouble(Column.ploidy)))
                    .duplicateProportion(Objects.requireNonNull(row.getDouble(Column.duplicateProportion)))
                    .meanReadsPerKb(Objects.requireNonNull(row.getDouble(Column.meanReadsPerKb)))
                    .gc50ReadsPerKb(Objects.requireNonNull(row.getDouble(Column.gc50ReadsPerKb)))
                .build()).collect(Collectors.toList());

            Validate.isTrue(telomereLengths.size() == 1, "Number of telomere length record must be 1");

            return telomereLengths.get(0);
        }
    }
}
