package com.hartwig.hmftools.tars.fasta;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;

// Concatenates a transcript's exons into a forward-strand sequence with its genomic exon spans.
public class TranscriptContigBuilder
{
    private final RefGenomeInterface mRefGenome;

    public TranscriptContigBuilder(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public TranscriptContigResult build(final GeneData gene, final TranscriptData transcript)
    {
        if(transcript.exons() == null || transcript.exons().isEmpty())
        {
            return null;
        }

        List<BaseRegion> spans = sortedExonSpans(transcript);

        StringBuilder sequence = new StringBuilder();
        int prevEnd = -1;
        for(BaseRegion span : spans)
        {
            if(span.start() <= prevEnd)
            {
                throw new IllegalStateException(String.format(
                        "overlapping exons in transcript %s gene %s: exon [%d,%d] overlaps prior exon ending at %d",
                        transcript.TransName, gene.GeneId, span.start(), span.end(), prevEnd));
            }

            sequence.append(mRefGenome.getBaseString(gene.Chromosome, span.start(), span.end()));
            prevEnd = span.end();
        }

        return new TranscriptContigResult(
                gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand,
                spans, sequence.toString());
    }

    // exons are stored by transcription rank; re-sort by genomic start for forward-strand orientation
    static List<BaseRegion> sortedExonSpans(final TranscriptData transcript)
    {
        List<ExonData> ordered = new ArrayList<>(transcript.exons());
        ordered.sort(Comparator.comparingInt(exon -> exon.Start));

        List<BaseRegion> spans = new ArrayList<>(ordered.size());
        for(ExonData exon : ordered)
        {
            spans.add(new BaseRegion(exon.Start, exon.End));
        }
        return spans;
    }

    public record TranscriptContigResult(
            String geneId,
            String geneName,
            String transName,
            String chromosome,
            int strand,
            List<BaseRegion> exonSpans,
            String sequence)
    {
    }
}
