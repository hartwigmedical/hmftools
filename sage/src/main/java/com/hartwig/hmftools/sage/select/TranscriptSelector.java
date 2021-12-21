package com.hartwig.hmftools.sage.select;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class TranscriptSelector
{
    private final List<TranscriptData> mTranscripts;
    private int mIndex = 0;

    public TranscriptSelector(final List<TranscriptData> regions)
    {
        mTranscripts = regions;
    }

    public TranscriptData select(int position)
    {
        if(mTranscripts.isEmpty())
            return null;

        TranscriptData current = mTranscripts.get(mIndex);

        while(mIndex > 0 && current.TransStart > position)
        {
            mIndex--;
            current = mTranscripts.get(mIndex);
        }

        while(mIndex < mTranscripts.size() - 1 && current.TransEnd < position)
        {
            mIndex++;
            current = mTranscripts.get(mIndex);
        }

        if(positionWithin(position, current.TransStart, current.TransEnd))
            return current;

        return null;
    }

}
