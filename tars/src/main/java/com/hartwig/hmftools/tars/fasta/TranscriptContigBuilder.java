package com.hartwig.hmftools.tars.fasta;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;

// Concatenates a transcript's exons into a forward-strand sequence with genomic span metadata.
// Placement onto the alt contig (altStart/altEnd) is SpliceFastaBuilder's responsibility.
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
            return null;

        // exons are stored by transcription rank; re-sort by genomic start for forward-strand orientation
        List<ExonData> ordered = new ArrayList<>(transcript.exons());
        ordered.sort(Comparator.comparingInt(e -> e.Start));

        List<BaseRegion> spans = new ArrayList<>(ordered.size());
        StringBuilder sequence = new StringBuilder();
        int prevEnd = -1;
        for(ExonData exon : ordered)
        {
            if(exon.Start <= prevEnd)
            {
                throw new IllegalStateException(String.format(
                        "overlapping exons in transcript %s gene %s: exon [%d,%d] overlaps prior exon ending at %d",
                        transcript.TransName, gene.GeneId, exon.Start, exon.End, prevEnd));
            }

            spans.add(new BaseRegion(exon.Start, exon.End));
            sequence.append(mRefGenome.getBaseString(gene.Chromosome, exon.Start, exon.End));
            prevEnd = exon.End;
        }

        return new TranscriptContigResult(
                gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, gene.Strand,
                spans, sequence.toString());
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
