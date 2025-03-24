package com.hartwig.hmftools.pavereverse.parse;

import static java.lang.String.format;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.TranscriptFilter;

class FilterTranscripts implements TranscriptRetriever
{
    private final ProteinVariantParser mVariantParser;
    private final boolean mReturnVariantForEachNonCanonicalTranscript;

    FilterTranscripts(final ProteinVariantParser proteinVariantParser, final boolean returnVariantForEachNonCanonicalTranscript)
    {
        mVariantParser = proteinVariantParser;
        mReturnVariantForEachNonCanonicalTranscript = returnVariantForEachNonCanonicalTranscript;
    }

    @Override
    public Set<ProteinTranscript> getApplicableTranscripts(GeneData geneData, TranscriptFilter refFilter)
    {
        List<TranscriptData> allTranscripts = mVariantParser.EnsemblCache.getTranscripts(geneData.GeneId);
        Set<ProteinTranscript> matchingRef = filter(allTranscripts, refFilter);
        if(matchingRef.isEmpty())
        {
            String msg = format("No transcript found for gene: %s matching ref: %s", geneData.GeneId, refFilter);
            throw new IllegalArgumentException(msg);
        }
        ProteinTranscript canonical = matchingRef.stream()
                .filter(pt -> pt.mTranscriptData.IsCanonical)
                .findFirst()
                .orElse(null);
        if(canonical != null)
        {
            return Set.of(canonical);
        }
        if(matchingRef.size() > 1 && !mReturnVariantForEachNonCanonicalTranscript)
        {
            String msg =
                    format("No canonical transcript, but multiple non-canonical transcripts, found for gene: %s, ref: %s", geneData.GeneId, refFilter);
            throw new IllegalArgumentException(msg);

        }
        return new HashSet<>(matchingRef);
    }

    private Set<ProteinTranscript> filter(List<TranscriptData> transcriptDataList, TranscriptFilter filter)
    {
        Set<ProteinTranscript> result = new HashSet<>();
        transcriptDataList.forEach(transcriptData ->
        {
            TranscriptAminoAcids aminoAcidsSequence = mVariantParser.TranscriptAminoAcidsMap.get(transcriptData.TransName);
            if(aminoAcidsSequence != null && filter.applies(aminoAcidsSequence))
            {
                result.add(new ProteinTranscript(aminoAcidsSequence, transcriptData));
            }
        });
        return result;
    }
}
