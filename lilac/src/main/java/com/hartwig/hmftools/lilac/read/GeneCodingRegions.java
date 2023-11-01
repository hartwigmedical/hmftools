package com.hartwig.hmftools.lilac.read;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;

public class GeneCodingRegions
{
    public final String GeneName;
    public final String Chromosome;
    public final byte Strand;
    public final int CodingStart;
    public final int CodingEnd;
    public final List<BaseRegion> CodingRegions;

    public GeneCodingRegions(final String geneName, final String chromosome, final TranscriptData transcriptData)
    {
        GeneName = geneName;
        Chromosome = chromosome;
        Strand = transcriptData.Strand;

        CodingStart = transcriptData.CodingStart;
        CodingEnd = transcriptData.CodingEnd;

        CodingRegions = Lists.newArrayList();

        int codingStart = transcriptData.CodingStart;
        int codingEnd = transcriptData.CodingEnd;

        for(ExonData exon : transcriptData.exons())
        {
            if(codingStart <= exon.End && codingEnd >= exon.Start)
            {
                CodingRegions.add(new BaseRegion(max(codingStart, exon.Start), min(codingEnd, exon.End)));
            }
        }
    }

    public boolean withinCodingBounds(int position, int buffer)
    {
        return position >= CodingStart - buffer && position <= CodingEnd + buffer;
    }

    public String toString()
    {
        return format("gene(%s) strand(%d) coding(%d-%d)", GeneName, Strand, CodingStart, CodingEnd);
    }
}
