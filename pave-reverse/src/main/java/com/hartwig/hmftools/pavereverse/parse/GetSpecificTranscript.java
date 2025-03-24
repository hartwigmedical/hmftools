package com.hartwig.hmftools.pavereverse.parse;

import static java.lang.String.format;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.TranscriptFilter;

class GetSpecificTranscript implements TranscriptRetriever
{
    private final ProteinVariantParser mVariantParser;
    private final String mTranscriptId;

    GetSpecificTranscript(final ProteinVariantParser proteinVariantParser, final String transcriptId)
    {
        mVariantParser = proteinVariantParser;
        mTranscriptId = transcriptId;
    }

    @Override
    public Set<ProteinTranscript> getApplicableTranscripts(final GeneData geneData, final TranscriptFilter refFilter)
    {
        TranscriptData transcriptData = mVariantParser.EnsemblCache.getTranscriptData(geneData.GeneId, mTranscriptId);
        if(transcriptData == null)
        {
            String msg = format("No transcript found. Gene: %s, transcript id: %s", geneData.GeneId, mTranscriptId);
            throw new IllegalArgumentException(msg);
        }
        TranscriptAminoAcids aminoAcidsSequence = mVariantParser.TranscriptAminoAcidsMap.get(transcriptData.TransName);
        if(!refFilter.applies(aminoAcidsSequence))
        {
            String m = format("Transcript does not match. Gene: %s, transcript: %s, ref: %s", geneData.GeneId, mTranscriptId, refFilter);
            throw new IllegalArgumentException(m);
        }
        return Set.of(new ProteinTranscript(aminoAcidsSequence, transcriptData));
    }
}
