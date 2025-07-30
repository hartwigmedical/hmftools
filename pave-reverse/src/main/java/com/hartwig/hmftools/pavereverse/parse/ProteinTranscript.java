package com.hartwig.hmftools.pavereverse.parse;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;

class ProteinTranscript
{
    final TranscriptAminoAcids mAminoAcids;
    final TranscriptData mTranscriptData;

    ProteinTranscript(TranscriptAminoAcids aminoAcids, TranscriptData transcriptData)
    {
        mAminoAcids = aminoAcids;
        mTranscriptData = transcriptData;
    }
}
