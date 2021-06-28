package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.Optional;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class TranscriptRegionSelector
{
    private final List<HmfTranscriptRegion> mRegions;
    private int mIndex = 0;

    public TranscriptRegionSelector(final List<HmfTranscriptRegion> regions)
    {
        mRegions = regions;
    }

    @NotNull
    public Optional<HmfTranscriptRegion> select(long position)
    {
        if(mRegions.isEmpty())
        {
            return Optional.empty();
        }

        HmfTranscriptRegion current = current();
        while(mIndex > 0 && current.start() > position)
        {
            mIndex--;
            current = current();
        }

        while(mIndex < mRegions.size() - 1 && current.end() < position)
        {
            mIndex++;
            current = current();
        }

        if(position >= current.start() && position <= current.end())
        {
            return Optional.of(current);
        }

        return Optional.empty();
    }

    private HmfTranscriptRegion current()
    {
        return mRegions.get(mIndex);
    }

}
