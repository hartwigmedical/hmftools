package com.hartwig.hmftools.geneutils.paneldesign;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;

public record GeneTranscript(
        GeneData gene,
        TranscriptData transcript
)
{
}
