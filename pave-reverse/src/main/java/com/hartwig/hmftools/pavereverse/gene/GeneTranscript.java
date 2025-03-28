package com.hartwig.hmftools.pavereverse.gene;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class GeneTranscript
{
    public final GeneData Gene;
    public final TranscriptData Transcript;
    public final List<ChrBaseRegion> CodingRegions;
    private final List<AnnotatedExon> mAnnotatedExons;
    private final int mIndexOfFirstTranslatedBase;
    private final int mIndexOfLastTranslatedBase;

    public GeneTranscript(final GeneData geneData, final TranscriptData transcriptData)
    {
        Gene = geneData;
        Transcript = transcriptData;
        mAnnotatedExons = new ArrayList<>();
        List<ChrBaseRegion> codingRegions = new ArrayList<>();
        final AtomicInteger basesCovered = new AtomicInteger(0);
        final AtomicBoolean firstCodingBaseFound = new AtomicBoolean(false);
        final AtomicInteger indexOfFirstCodingBase = new AtomicInteger(0);
        final AtomicInteger indexOfLastCodingBase = new AtomicInteger(0);
        if(transcriptData.posStrand())
        {
            transcriptData.exons().forEach(exon ->
            {
                int firstBaseInExon = basesCovered.get() + 1;
                int lastBaseInExon = firstBaseInExon + exon.baseLength() - 1;
                mAnnotatedExons.add(new AnnotatedExon(firstBaseInExon, lastBaseInExon, exon));
                if(exon.End >= transcriptData.CodingStart && exon.Start <= transcriptData.CodingEnd)
                {
                    if(!firstCodingBaseFound.get())
                    {
                        indexOfFirstCodingBase.set(firstBaseInExon + (transcriptData.CodingStart - exon.Start));
                        firstCodingBaseFound.set(true);
                    }
                    if(exon.End >= transcriptData.CodingEnd)
                    {
                        indexOfLastCodingBase.set(firstBaseInExon + (transcriptData.CodingEnd - exon.Start));
                    }
                    codingRegions.add(new ChrBaseRegion(geneData.Chromosome, Math.max(exon.Start, transcriptData.CodingStart), Math.min(exon.End, transcriptData.CodingEnd)));
                }
                basesCovered.set(lastBaseInExon);
            });
        }
        else
        {
            List<ExonData> exons = new ArrayList<>(transcriptData.exons());
            Collections.reverse(exons);
            exons.forEach(exon ->
            {
                int firstBaseInExon = basesCovered.get() + 1;
                int lastBaseInExon = firstBaseInExon + exon.baseLength() - 1;
                mAnnotatedExons.add(new AnnotatedExon(firstBaseInExon, lastBaseInExon, exon, false));
                if(exon.End >= transcriptData.CodingStart && exon.Start <= transcriptData.CodingEnd)
                {
                    if(!firstCodingBaseFound.get())
                    {
                        indexOfFirstCodingBase.set(firstBaseInExon + (exon.End - transcriptData.CodingEnd));
                        firstCodingBaseFound.set(true);
                    }
                    if(exon.Start <= transcriptData.CodingStart)
                    {
                        indexOfLastCodingBase.set(firstBaseInExon + (exon.End - transcriptData.CodingStart));
                    }
                    codingRegions.add(new ChrBaseRegion(geneData.Chromosome, Math.max(exon.Start, transcriptData.CodingStart), Math.min(exon.End, transcriptData.CodingEnd)));
                }
                basesCovered.set(lastBaseInExon);
            });
        }
        CodingRegions = codingRegions;
        mIndexOfFirstTranslatedBase = indexOfFirstCodingBase.get();
        mIndexOfLastTranslatedBase = indexOfLastCodingBase.get();
        Preconditions.checkState(totalTranslatedLength() == (mIndexOfLastTranslatedBase - mIndexOfFirstTranslatedBase + 1));
    }

    public List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
    }

    public boolean isCodingBaseAtTheStartOfAnExon(int base)
    {
        int baseToSeek = base + mIndexOfFirstTranslatedBase - 1;
        return mAnnotatedExons.stream().anyMatch(annotatedExon -> annotatedExon.FirstBase == baseToSeek);
    }

    public boolean isCodingBaseAtTheEndOfAnExon(int base)
    {
        int baseToSeek = base + mIndexOfFirstTranslatedBase - 1;
        return mAnnotatedExons.stream().anyMatch(annotatedExon -> annotatedExon.LastBase == baseToSeek);
    }

    public int totalTranslatedLength()
    {
        return CodingRegions.stream().mapToInt(ChrBaseRegion::baseLength).sum();
    }

    public int absolutePositionOfTranslatedBase(int baseIndex)
    {
        Preconditions.checkArgument(baseIndex >= 1);
        Preconditions.checkArgument(baseIndex <= totalTranslatedLength());
        return absolutePosition(baseIndex + mIndexOfFirstTranslatedBase - 1);
    }

    public int absolutePositionOf5PrimeUtrExonicBase(int baseIndex)
    {
        Preconditions.checkArgument(baseIndex < 0);
        return absolutePosition(baseIndex + mIndexOfFirstTranslatedBase);
    }

    public int absolutePositionOf3PrimeUtrExonicBase(int baseIndex)
    {
        Preconditions.checkArgument(baseIndex > 0);
        return absolutePosition(baseIndex + mIndexOfLastTranslatedBase);
    }

    private AnnotatedExon seek(int baseIndex)
    {
        return mAnnotatedExons.stream().filter(annotatedExon -> annotatedExon.contains(baseIndex)).findFirst().orElse(null);
    }

    private int absolutePosition(int index)
    {
        if(index <= 0)
        {
            return mAnnotatedExons.get(0).getAbsolutePositionOfBaseUpstreamOfExon(index);
        }
        AnnotatedExon exon = seek(index);
        if(exon == null)
        {
            throw new IllegalArgumentException("No exon found for index " + index);
        }
        return exon.getAbsolutePosition(index);
    }
}
