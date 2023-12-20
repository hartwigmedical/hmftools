package com.hartwig.hmftools.common.purple;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public final class HrdDataFile
{
    private static final String EXTENSION = ".purple.hrd.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return basePath + File.separator + sample + EXTENSION;
    }

    enum Column
    {
        lohSegments,
        segmentImbalances,
        segmentBreaks,
        status;
    }

    public static HrdData read(final String basePath, final String sample) throws IOException
    {
        final String filename = generateFilename(basePath, sample);

        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            for(DelimFileReader.Row row : reader)
            {
                int lohSegments = row.getInt(Column.lohSegments);
                int segmentImbalances = row.getInt(Column.segmentImbalances);
                int segmentBreaks = row.getInt(Column.segmentBreaks);

                HrdStatus status = HrdStatus.valueOf(row.get(Column.status));
                return new HrdData(lohSegments, segmentImbalances, segmentBreaks, status);
            }
        }

        return null;
    }

    public static void write(final String basePath, final String sample, final HrdData hrdData) throws IOException
    {
        final String fileName = generateFilename(basePath, sample);

        try(BufferedWriter writer = createBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, Column.values(), List.of(hrdData),
                    (hrd, row) -> {
                        row.set(Column.lohSegments, hrd.LohSegments);
                        row.set(Column.segmentImbalances, hrd.SegmentImbalances);
                        row.set(Column.segmentBreaks, hrd.SegmentBreaks);
                        row.set(Column.status, hrd.Status.toString());
                    });
        }
    }
}
