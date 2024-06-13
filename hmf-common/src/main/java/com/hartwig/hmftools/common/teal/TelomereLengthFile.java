package com.hartwig.hmftools.common.teal;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.Validate;

public class TelomereLengthFile
{
    enum Column
    {
        sampleId, type, rawTelomereLength, finalTelomereLength, fullFragments, cRichPartialFragments, gRichPartialFragments,
        totalTelomericReads, purity, ploidy, duplicateProportion, meanReadDepth, gc50ReadDepth
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
                    .type(row.get(Column.type))
                    .rawTelomereLength(row.getDouble(Column.rawTelomereLength))
                    .finalTelomereLength(row.getDouble(Column.finalTelomereLength))
                    .fullFragments(row.getInt(Column.fullFragments))
                    .cRichPartialFragments(row.getInt(Column.cRichPartialFragments))
                    .gRichPartialFragments(row.getInt(Column.gRichPartialFragments))
                    .totalTelomericReads(row.getInt(Column.totalTelomericReads))
                    .purity(row.getDouble(Column.purity))
                    .ploidy(row.getDouble(Column.ploidy))
                    .duplicateProportion(row.getDouble(Column.duplicateProportion))
                    .meanReadDepth(row.getDouble(Column.meanReadDepth))
                    .gc50ReadDepth(row.getDouble(Column.gc50ReadDepth))
                .build()).collect(Collectors.toList());

            Validate.isTrue(telomereLengths.size() == 1, "Number of telomere length record must be 1");

            return telomereLengths.get(0);
        }
    }

    public static void write(final String filename, final String sampleId, TelomereLength telomereLength)
    {
        try (DelimFileWriter<TelomereLength> writer = new DelimFileWriter<>(filename, Column.values(),
                (tl, row) ->
                {
                    row.set(Column.sampleId, sampleId);
                    row.set(Column.type, tl.type());
                    row.set(Column.rawTelomereLength, tl.rawTelomereLength(), "%.2f");
                    row.set(Column.finalTelomereLength, tl.finalTelomereLength(), "%.2f");
                    row.set(Column.fullFragments, tl.fullFragments());
                    row.set(Column.cRichPartialFragments, tl.cRichPartialFragments());
                    row.set(Column.gRichPartialFragments, tl.gRichPartialFragments());
                    row.set(Column.totalTelomericReads, tl.totalTelomericReads());
                    row.set(Column.purity, tl.purity());
                    row.set(Column.ploidy, tl.ploidy());
                    row.set(Column.duplicateProportion, tl.duplicateProportion());
                    row.set(Column.meanReadDepth, tl.meanReadDepth());
                    row.set(Column.gc50ReadDepth, tl.gc50ReadDepth());
                }))
        {
            writer.writeRow(telomereLength);
        }
    }
}
