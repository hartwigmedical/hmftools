package com.hartwig.hmftools.pavereverse.gene;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class GeneTranscript
{
    public final GeneData Gene;
    public final TranscriptData Transcript;
    public final List<ChrBaseRegion> CodingRegions;

    public GeneTranscript(final GeneData geneData, final TranscriptData transcriptData)
    {
        this.Gene = geneData;
        this.Transcript = transcriptData;
        List<ChrBaseRegion> codingRegions = transcriptData.exons().stream()
                .filter(exonData -> exonData.End >= transcriptData.CodingStart && exonData.Start <= transcriptData.CodingEnd)
                .map(exonData -> new ChrBaseRegion(geneData.Chromosome, Math.max(exonData.Start, transcriptData.CodingStart), Math.min(exonData.End, transcriptData.CodingEnd)))
                .collect(Collectors.toList());
        if(transcriptData.negStrand())
        {
            List<ChrBaseRegion> reversed = new ArrayList<>(codingRegions);
            Collections.reverse(reversed);
            CodingRegions = Collections.unmodifiableList(reversed);
        }
        else
        {
            CodingRegions = codingRegions;
        }
    }

    public List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
    }
}
