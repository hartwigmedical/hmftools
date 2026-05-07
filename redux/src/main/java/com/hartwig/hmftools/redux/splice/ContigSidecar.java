package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public final class ContigSidecar
{
    public static final String COL_CONTIG_NAME = "ContigName";
    public static final String COL_GENE_ID = "GeneId";
    public static final String COL_GENE_NAME = "GeneName";
    public static final String COL_TRANS_NAME = "TransName";
    public static final String COL_CHROMOSOME = "Chromosome";
    public static final String COL_EXON_SPANS = "ExonSpans";

    private static final List<String> COLUMNS = List.of(
            COL_CONTIG_NAME, COL_GENE_ID, COL_GENE_NAME, COL_TRANS_NAME, COL_CHROMOSOME, COL_EXON_SPANS);

    // matches the convention in SpecificRegions: ITEM_DELIM separates spans, '-' separates start from end (BaseRegion.toString)
    private static final String SPAN_DELIM = "-";

    public static void write(final String filename, final List<ContigEntry> entries)
    {
        try(DelimFileWriter<ContigEntry> writer = new DelimFileWriter<>(filename, COLUMNS, (entry, row) ->
        {
            row.set(COL_CONTIG_NAME, entry.ContigName);
            row.set(COL_GENE_ID, entry.GeneId);
            row.set(COL_GENE_NAME, entry.GeneName);
            row.set(COL_TRANS_NAME, entry.TransName);
            row.set(COL_CHROMOSOME, entry.Chromosome);
            row.set(COL_EXON_SPANS, entry.ExonSpans.stream().map(BaseRegion::toString).collect(Collectors.joining(ITEM_DELIM)));
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
                        row.get(COL_CONTIG_NAME),
                        row.get(COL_GENE_ID),
                        row.get(COL_GENE_NAME),
                        row.get(COL_TRANS_NAME),
                        row.get(COL_CHROMOSOME),
                        decodeSpans(row.get(COL_EXON_SPANS))));
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
