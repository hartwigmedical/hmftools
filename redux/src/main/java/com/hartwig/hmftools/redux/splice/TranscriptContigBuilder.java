package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.splice.SpliceCommon.CONTIG_NAME_DELIM;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.CONTIG_NAME_PREFIX;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.BaseRegion;

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
        for(ExonData exon : ordered)
        {
            spans.add(new BaseRegion(exon.Start, exon.End));
            sequence.append(mRefGenome.getBaseString(gene.Chromosome, exon.Start, exon.End));
        }

        String contigName = CONTIG_NAME_PREFIX + gene.GeneId
                + CONTIG_NAME_DELIM + gene.GeneName
                + CONTIG_NAME_DELIM + transcript.TransName;

        ContigEntry entry = new ContigEntry(
                contigName, gene.GeneId, gene.GeneName, transcript.TransName, gene.Chromosome, spans);

        return new TranscriptContigResult(entry, sequence.toString());
    }

    public static final class TranscriptContigResult
    {
        public final ContigEntry Entry;
        public final String Sequence;

        public TranscriptContigResult(final ContigEntry entry, final String sequence)
        {
            Entry = entry;
            Sequence = sequence;
        }
    }
}
