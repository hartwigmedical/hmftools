package com.hartwig.hmftools.pavereverse.gene;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

import com.google.common.base.Preconditions;
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
        final AtomicBoolean startFound = new AtomicBoolean(false);
        final AtomicInteger indexOfFirstTranslatedBase = new AtomicInteger(0);
        final AtomicInteger indexOfLastTranslatedBase = new AtomicInteger(0);
        transcriptData.exons().forEach(exon ->
        {
            int firstBaseInExon = basesCovered.get() + 1;
            int lastBaseInExon = firstBaseInExon + exon.baseLength() - 1;
            mAnnotatedExons.add(new AnnotatedExon(firstBaseInExon, lastBaseInExon, exon));
            if(exon.End >= transcriptData.CodingStart && exon.Start <= transcriptData.CodingEnd)
            {
                if(!startFound.get())
                {
                    indexOfFirstTranslatedBase.set(firstBaseInExon + (transcriptData.CodingStart - exon.Start));
                    startFound.set(true);
                }
                if(exon.End >= transcriptData.CodingEnd)
                {
                    indexOfLastTranslatedBase.set(firstBaseInExon + (transcriptData.CodingEnd - exon.Start));
                }
                codingRegions.add(new ChrBaseRegion(geneData.Chromosome, Math.max(exon.Start, transcriptData.CodingStart), Math.min(exon.End, transcriptData.CodingEnd)));
            }
            basesCovered.set(lastBaseInExon);
        });
        if(transcriptData.negStrand())
        {
            List<ChrBaseRegion> reversed = new ArrayList<>(codingRegions);
            Collections.reverse(reversed);
            CodingRegions = Collections.unmodifiableList(reversed);
            mIndexOfFirstTranslatedBase = -99; // todo
            mIndexOfLastTranslatedBase = -99;
        }
        else
        {
            CodingRegions = codingRegions;
            mIndexOfFirstTranslatedBase = indexOfFirstTranslatedBase.get();
            mIndexOfLastTranslatedBase = indexOfLastTranslatedBase.get();
            Preconditions.checkState(totalTranslatedLength() == (mIndexOfLastTranslatedBase - mIndexOfFirstTranslatedBase + 1));
        }
    }

    public List<Integer> codingRegionLengths()
    {
        return CodingRegions.stream().map(ChrBaseRegion::baseLength).collect(Collectors.toList());
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
        AnnotatedExon exon = seek(index);
        if(exon == null)
        {
            throw new IllegalArgumentException("No exon found for index " + index);
        }
        return exon.getAbsolutePosition(index);
    }
}
