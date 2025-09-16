package com.hartwig.hmftools.lilac.read;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.lilac.hla.HlaGene_;

public class GeneCodingRegions
{
    public final HlaGene_ GeneName;
    public final String Chromosome;
    public final byte Strand;
    public final int CodingStart;
    public final int CodingEnd;
    public final List<BaseRegion> CodingRegions;

    public GeneCodingRegions(final HlaGene_ geneName, final String chromosome, final TranscriptData transcriptData_)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Strand = transcriptData_.Strand;

        CodingStart = transcriptData_.CodingStart;
        CodingEnd = transcriptData_.CodingEnd;

        CodingRegions = Lists.newArrayList();

        int codingStart = transcriptData_.CodingStart;
        int codingEnd = transcriptData_.CodingEnd;

        for(ExonData exon_ : transcriptData_.exons())
        {
            if(codingStart <= exon_.End && codingEnd >= exon_.Start)
            {
                CodingRegions.add(new BaseRegion(max(codingStart, exon_.Start), min(codingEnd, exon_.End)));
            }
        }
    }

    public boolean withinCodingBounds(int position, int buffer)
    {
        return position >= CodingStart - buffer && position <= CodingEnd + buffer;
    }

    @Override
    public String toString()
    {
        return format("gene(%s) strand(%d) coding(%d-%d)", GeneName, Strand, CodingStart, CodingEnd);
    }
}
