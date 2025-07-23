package com.hartwig.hmftools.pavereverse.parse;

import java.util.Set;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.pavereverse.TranscriptFilter;

interface TranscriptRetriever
{
    Set<ProteinTranscript> getApplicableTranscripts(GeneData geneData, TranscriptFilter refFilter);
}
