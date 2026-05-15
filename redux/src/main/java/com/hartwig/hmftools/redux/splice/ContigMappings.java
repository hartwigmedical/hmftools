package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public final class ContigMappings
{
    public enum Column
    {
        ContigName,
        GeneId,
        GeneName,
        TransName,
        Chromosome,
        ExonSpans
    }

    // matches the convention in SpecificRegions: ITEM_DELIM separates spans, '-' separates start from end (BaseRegion.toString)
    private static final String SPAN_DELIM = "-";

    public static void write(final String filename, final List<ContigEntry> entries)
    {
        try(DelimFileWriter<ContigEntry> writer = new DelimFileWriter<>(filename, Column.values(), (entry, row) ->
        {
            row.set(Column.ContigName, entry.contigName());
            row.set(Column.GeneId, entry.geneId());
            row.set(Column.GeneName, entry.geneName());
            row.set(Column.TransName, entry.transName());
            row.set(Column.Chromosome, entry.chromosome());
            row.set(Column.ExonSpans, entry.exonSpans().stream().map(BaseRegion::toString).collect(Collectors.joining(ITEM_DELIM)));
        }))
        {
            for(ContigEntry entry : entries)
                writer.writeRow(entry);
        }

        RD_LOGGER.info("wrote {} contig entries to {}", entries.size(), filename);
    }

    public static List<ContigEntry> read(final String filename)
    {
        List<ContigEntry> entries = new ArrayList<>();

        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            for(DelimFileReader.Row row : reader)
            {
                entries.add(new ContigEntry(
                        row.get(Column.ContigName),
                        row.get(Column.GeneId),
                        row.get(Column.GeneName),
                        row.get(Column.TransName),
                        row.get(Column.Chromosome),
                        decodeSpans(row.get(Column.ExonSpans))));
            }
        }

        RD_LOGGER.info("loaded {} contig entries from {}", entries.size(), filename);
        return entries;
    }

    private static List<BaseRegion> decodeSpans(final String encoded)
    {
        List<BaseRegion> spans = new ArrayList<>();
        for(String token : encoded.split(ITEM_DELIM))
        {
            String[] parts = token.split(SPAN_DELIM);
            spans.add(new BaseRegion(Integer.parseInt(parts[0]), Integer.parseInt(parts[1])));
        }
        return spans;
    }

}
