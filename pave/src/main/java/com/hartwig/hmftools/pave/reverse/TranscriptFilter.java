package com.hartwig.hmftools.pave.reverse;

import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;

interface TranscriptFilter
{
    boolean applies(TranscriptAminoAcids aminoAcids);
}
