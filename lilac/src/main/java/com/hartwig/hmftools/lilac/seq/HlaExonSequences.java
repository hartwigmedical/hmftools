package com.hartwig.hmftools.lilac.seq;

import static java.lang.Math.ceil;
import static java.lang.Math.min;

import static com.hartwig.hmftools.lilac.LilacConstants.EXON_CHUNK_SIZE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;

public final class HlaExonSequences
{
    public record ExonSequence(List<String> sequences, boolean hasWildcards) {}

    public final HlaAllele Allele;
    public final ImmutableList<ExonSequence> ExonSequences;

    private HlaExonSequences(final HlaAllele allele, final Iterable<ExonSequence> exonSequences)
    {
        Allele = allele;
        ExonSequences = ImmutableList.copyOf(exonSequences);
    }

    public static HlaExonSequences create(final Map<HlaGene, List<Integer>> geneExonBoundaries, final HlaSequenceLoci sequence)
    {
        HlaAllele allele = sequence.Allele;
        List<Integer> exonBoundaries = geneExonBoundaries.get(allele.Gene);
        List<String> acids = sequence.getSequences();
        List<List<String>> exonAcids = Lists.newArrayList();

        int index = 0;
        for(int exonBoundary : exonBoundaries)
        {
            if(index >= acids.size())
                break;

            int toIndex = min(exonBoundary + 1, acids.size());
            List<String> exonAcid = acids.subList(index, toIndex);

            List<List<String>> chunks = Lists.newArrayList();
            int chunkCount = (int) ceil(((double) exonAcid.size()) / EXON_CHUNK_SIZE);
            int chunkSize = (int) ceil(((double) exonAcid.size()) / chunkCount);
            int start = 0;
            while(start < exonAcid.size())
            {
                int end = min(start + chunkSize, exonAcid.size());
                chunks.add(exonAcid.subList(start, end));
                start = end;
            }

            exonAcids.addAll(chunks);
            index = exonBoundary + 1;
        }

        List<ExonSequence> exonSeqs = Lists.newArrayList();
        for(List<String> seq : exonAcids)
        {
            boolean hasWildcards = seq.stream().anyMatch("*"::equals);
            exonSeqs.add(new ExonSequence(seq, hasWildcards));
        }

        return new HlaExonSequences(allele, exonSeqs);
    }
}
