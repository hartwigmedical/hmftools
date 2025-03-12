package com.hartwig.hmftools.pavereverse;

import java.util.Set;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

public final class BaseSequenceVariants
{
    public final TranscriptData Transcript;
    public final String Chromosome;
    final Set<BaseSequenceChange> Changes;

    public BaseSequenceVariants(TranscriptData transcript, String chromosome, Set<BaseSequenceChange> changes)
    {
        Transcript = transcript;
        Chromosome = RefGenomeFunctions.stripChrPrefix(chromosome);
        Changes = changes;
    }

    public String transcriptName()
    {
        return Transcript.TransName;
    }

    public Set<BaseSequenceChange> changes()
    {
        return Changes;
    }
}
