package com.hartwig.hmftools.pavereverse;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

public interface TranscriptFilter
{
    boolean applies(TranscriptAminoAcids aminoAcids);
}
