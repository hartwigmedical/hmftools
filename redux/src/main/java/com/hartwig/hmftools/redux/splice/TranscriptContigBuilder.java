package com.hartwig.hmftools.redux.splice;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;

// per-transcript exon-concatenation. Produces a forward-strand sequence and the genomic exon spans that
// generated it. Placement of this sequence onto a per-chromosome alt contig (altStart/altEnd) is the
// caller's job — SpliceFastaBuilder owns that packing.
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

        // exons are stored by Rank (transcribed order). Re-sort by genomic Start so the contig sequence is in
        // forward-strand orientation regardless of gene strand — keeps lift-back math direct.
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
                gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, spans, sequence.toString());
    }

    public static final class TranscriptContigResult
    {
        public final String GeneId;
        public final String GeneName;
        public final String TransName;
        public final String Chromosome;
        public final List<BaseRegion> ExonSpans;
        public final String Sequence;

        public TranscriptContigResult(
                final String geneId, final String geneName, final String transName, final String chromosome,
                final List<BaseRegion> exonSpans, final String sequence)
        {
            GeneId = geneId;
            GeneName = geneName;
            TransName = transName;
            Chromosome = chromosome;
            ExonSpans = exonSpans;
            Sequence = sequence;
        }
    }
}
