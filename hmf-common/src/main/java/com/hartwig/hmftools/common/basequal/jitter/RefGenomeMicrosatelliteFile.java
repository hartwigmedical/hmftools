package com.hartwig.hmftools.common.basequal.jitter;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.util.StringUtil;

public class RefGenomeMicrosatelliteFile implements AutoCloseable
{
    enum Column
    {
        chromosome,
        start,
        end,
        unit,
        mappability
    }

    private static final String FILE_EXTENSION = ".ref_genome_ms.tsv.gz";

    public static String generateFilename(String basePath, RefGenomeVersion refGenomeVersion)
    {
        return basePath + File.separator + refGenomeVersion + FILE_EXTENSION;
    }

    private final DelimFileWriter<RefGenomeMicrosatellite> mWriter;

    public RefGenomeMicrosatelliteFile(String fileName)
    {
        mWriter = new DelimFileWriter<>(fileName, Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList()),
                (refGenomeMicrosatellite, row) ->
                {
                    row.set(Column.chromosome, refGenomeMicrosatellite.chromosome());
                    row.set(Column.start, refGenomeMicrosatellite.referenceStart());
                    row.set(Column.end, refGenomeMicrosatellite.referenceEnd());
                    row.set(Column.unit, refGenomeMicrosatellite.unitString());
                    row.set(Column.mappability, refGenomeMicrosatellite.mappability);
                });
    }

    @Override
    public void close() throws Exception
    {
        mWriter.close();
    }

    public void writeRow(@NotNull final RefGenomeMicrosatellite refGenomeMicrosatellite)
    {
        mWriter.writeRow(refGenomeMicrosatellite);
    }

    public static List<RefGenomeMicrosatellite> read(final String filename)
    {
        List<RefGenomeMicrosatellite> refGenomeMicrosatellites = new ArrayList<>(100_000);
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            reader.stream().map(row ->
                        new RefGenomeMicrosatellite(row.get(Column.chromosome), row.getInt(Column.start), row.getInt(Column.end),
                                StringUtil.stringToBytes(row.get(Column.unit)))).forEach(refGenomeMicrosatellites::add);
        }
        return refGenomeMicrosatellites;
    }
}
