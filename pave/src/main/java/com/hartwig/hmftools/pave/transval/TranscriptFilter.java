package com.hartwig.hmftools.pave.transval;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

public interface TranscriptFilter
{
    boolean applies(TranscriptAminoAcids aminoAcids);
}
