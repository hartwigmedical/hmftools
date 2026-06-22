package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public final class ContigSidecar
{
    public enum Column
    {
        ContigName,
        AltStart,
        AltEnd,
        GeneId,
        GeneName,
        TransName,
        Chromosome,
        Strand,
        ExonSpans
    }

    private static final String SPAN_DELIM = "-";

    public static void write(final String filename, final List<ContigEntry> entries)
    {
        try(DelimFileWriter<ContigEntry> writer = new DelimFileWriter<>(filename, Column.values(), (entry, row) ->
        {
            row.set(Column.ContigName, entry.contigName());
            row.set(Column.AltStart, String.valueOf(entry.altStart()));
            row.set(Column.AltEnd, String.valueOf(entry.altEnd()));
            row.set(Column.GeneId, entry.geneId());
            row.set(Column.GeneName, entry.geneName());
            row.set(Column.TransName, entry.transName());
            row.set(Column.Chromosome, entry.chromosome());
            row.set(Column.Strand, String.valueOf(entry.strand()));
            row.set(Column.ExonSpans, entry.exonSpans().stream().map(BaseRegion::toString).collect(Collectors.joining(ITEM_DELIM)));
        }))
        {
            for(ContigEntry entry : entries)
            {
                writer.writeRow(entry);
            }
        }

        TARS_LOGGER.info("wrote {} contig entries to {}", entries.size(), filename);
    }

    public static List<ContigEntry> read(final String filename)
    {
        List<ContigEntry> entries = new ArrayList<>();
        boolean missingStrandWarned = false;

        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            // Strand was added after the initial release; older sidecars lack it. Default to 0 (suppresses XS:A emission) rather than failing.
            boolean hasStrand = reader.getColumnNames().contains(Column.Strand.name());

            for(DelimFileReader.Row row : reader)
            {
                int strand;
                if(hasStrand)
                {
                    strand = Integer.parseInt(row.get(Column.Strand));
                }
                else
                {
                    strand = 0;
                    if(!missingStrandWarned)
                    {
                        TARS_LOGGER.warn("contig sidecar {} lacks the Strand column - XS:A:+/- will not "
                                + "be emitted on tx-derived spliced records. Regenerate the sidecar "
                                + "with the current SpliceFastaBuilder to enable.", filename);
                        missingStrandWarned = true;
                    }
                }

                entries.add(new ContigEntry(
                        row.get(Column.ContigName),
                        Integer.parseInt(row.get(Column.AltStart)),
                        Integer.parseInt(row.get(Column.AltEnd)),
                        row.get(Column.GeneId),
                        row.get(Column.GeneName),
                        row.get(Column.TransName),
                        row.get(Column.Chromosome),
                        strand,
                        decodeSpans(row.get(Column.ExonSpans))));
            }
        }

        TARS_LOGGER.info("loaded {} contig entries from {}", entries.size(), filename);
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
