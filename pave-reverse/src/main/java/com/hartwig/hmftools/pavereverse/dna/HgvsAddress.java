package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public interface HgvsAddress
{
    int toStrandLocation(GeneTranscript geneTranscript);
}
