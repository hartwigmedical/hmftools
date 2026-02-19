package com.hartwig.hmftools.redux.jitter;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

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
        return format("%s/msi_jitter_sites.%s.tsv.gz", basePath, refGenomeVersion.identifier());
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

    public void writeRow(final MicrosatelliteSite microsatelliteSite)
    {
        mWriter.writeRow(microsatelliteSite);
    }

    public static List<MicrosatelliteSite> read(final String filename)
    {
        List<MicrosatelliteSite> microsatelliteSites = new ArrayList<>(3_000_000);

        try
        {
            BufferedReader reader = createBufferedReader(filename);

            String header = reader.readLine();
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(Column.chromosome.toString());
            int startIndex = fieldsIndexMap.get(Column.start.toString());
            int endIndex = fieldsIndexMap.get(Column.end.toString());
            int unitIndex = fieldsIndexMap.get(Column.unit.toString());

            String line = null;

            while((line = reader.readLine()) != null)
            {
                String[] values = line.split(TSV_DELIM, -1);

                MicrosatelliteSite msSite = new MicrosatelliteSite(
                        values[chrIndex],
                        Integer.parseInt(values[startIndex]),
                        Integer.parseInt(values[endIndex]),
                        StringUtil.stringToBytes(values[unitIndex]));

                microsatelliteSites.add(msSite);
            }

            return microsatelliteSites;
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to load MS indel sites file({}): {}", filename, e.toString());
            return null;
        }
    }
}
