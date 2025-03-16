package com.hartwig.hmftools.pavereverse.parse;

import static java.lang.String.format;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.TranscriptFilter;

class GetSpecificTranscript implements TranscriptRetriever
{
    private final ProteinVariantParser proteinVariantParser;
    final String transcriptId;

    GetSpecificTranscript(final ProteinVariantParser proteinVariantParser, final String transcriptId)
    {
        this.proteinVariantParser = proteinVariantParser;
        this.transcriptId = transcriptId;
    }

    @Override
    public Set<ProteinTranscript> getApplicableTranscripts(final GeneData geneData, final TranscriptFilter refFilter)
    {
        TranscriptData transcriptData = proteinVariantParser.EnsemblCache.getTranscriptData(geneData.GeneId, transcriptId);
        if(transcriptData == null)
        {
            String msg = format("No transcript found. Gene: %s, transcript id: %s", geneData.GeneId, transcriptId);
            throw new IllegalArgumentException(msg);
        }
        TranscriptAminoAcids aminoAcidsSequence = proteinVariantParser.TranscriptAminoAcidsMap.get(transcriptData.TransName);
        if(!refFilter.applies(aminoAcidsSequence))
        {
            String m = format("Transcript does not match. Gene: %s, transcript: %s, ref: %s", geneData.GeneId, transcriptId, refFilter);
            throw new IllegalArgumentException(m);
        }
        return Set.of(new ProteinTranscript(aminoAcidsSequence, transcriptData));
    }
}
