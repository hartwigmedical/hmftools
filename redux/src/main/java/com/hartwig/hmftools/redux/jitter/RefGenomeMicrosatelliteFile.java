package com.hartwig.hmftools.redux.jitter;

import static java.lang.String.format;

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

    public static String generateFilename(final String basePath, final RefGenomeVersion refGenomeVersion)
    {
        return format("%s/msi_jitter_sites.%d.tsv.gz", basePath, refGenomeVersion.identifier());
    }

    private final DelimFileWriter<MicrosatelliteSite> mWriter;

    public RefGenomeMicrosatelliteFile(String fileName)
    {
        mWriter = new DelimFileWriter<>(fileName, Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList()),
                (refGenomeMicrosatellite, row) ->
                {
                    row.set(Column.chromosome, refGenomeMicrosatellite.chromosome());
                    row.set(Column.start, refGenomeMicrosatellite.referenceStart());
                    row.set(Column.end, refGenomeMicrosatellite.referenceEnd());
                    row.set(Column.unit, refGenomeMicrosatellite.unitString());
                    row.set(Column.mappability, refGenomeMicrosatellite.mappability());
                });
    }

    @Override
    public void close()
    {
        mWriter.close();
    }

    public void writeRow(@NotNull final MicrosatelliteSite microsatelliteSite)
    {
        mWriter.writeRow(microsatelliteSite);
    }

    public static List<MicrosatelliteSite> read(final String filename)
    {
        List<MicrosatelliteSite> microsatelliteSites = new ArrayList<>(100_000);
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            reader.stream().map(row ->
                        new MicrosatelliteSite(
                                row.get(Column.chromosome),
                                row.getInt(Column.start),
                                row.getInt(Column.end),
                                StringUtil.stringToBytes(row.get(Column.unit)))).forEach(microsatelliteSites::add);
        }
        return microsatelliteSites;
    }
}
