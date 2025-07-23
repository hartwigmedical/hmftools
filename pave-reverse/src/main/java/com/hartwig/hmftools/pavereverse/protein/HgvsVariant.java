package com.hartwig.hmftools.pavereverse.protein;

import java.util.List;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public abstract class HgvsVariant
{
    private final GeneTranscript mGeneTranscript;

    protected HgvsVariant(GeneData gene, TranscriptData transcript)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        mGeneTranscript = new GeneTranscript(gene, transcript);
    }

    public String geneName()
    {
        return mGeneTranscript.Gene.GeneName;
    }

    public String geneId()
    {
        return mGeneTranscript.Gene.GeneId;
    }

    public String transcriptName()
    {
        return mGeneTranscript.Transcript.TransName;
    }

    public boolean forwardStrand()
    {
        return mGeneTranscript.Gene.forwardStrand();
    }

    public boolean reverseStrand()
    {
        return mGeneTranscript.Gene.reverseStrand();
    }

    public String chromosome()
    {
        return mGeneTranscript.Gene.Chromosome;
    }

    public List<ChrBaseRegion> codingRegions()
    {
        return mGeneTranscript.CodingRegions;
    }

    public TranscriptData transcriptData()
    {
        return mGeneTranscript.Transcript;
    }

    public GeneTranscript geneTranscript()
    {
        return mGeneTranscript;
    }

    @VisibleForTesting
    List<Integer> codingRegionLengths()
    {
        return mGeneTranscript.codingRegionLengths();
    }
}
