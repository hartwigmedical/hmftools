package com.hartwig.hmftools.svprep.reads;

public class ReadIdTrimmer
{
    private final boolean mEnabled;
    private int mReadIdTrimIndex;

    public ReadIdTrimmer(boolean enabled)
    {
        mEnabled = enabled;
        mReadIdTrimIndex = -1;
    }

    public String trim(final String readId)
    {
        if(!mEnabled)
            return readId;

        if(mReadIdTrimIndex < 0)
            mReadIdTrimIndex = findReadIdTrimIndex(readId);

        return mReadIdTrimIndex > 0 ? readId.substring(mReadIdTrimIndex) : readId;
    }

    private static final char READ_ID_DELIM = ':';

    private static int findReadIdTrimIndex(final String readId)
    {
        // expected IDs: A00260:251:HLYGFDSXY:1:1673:32280:4946, ultima 011852_2-UGAv3-2-1333458495
        int delimCount = 0;
        for(int i = 0; i < readId.length(); ++i)
        {
            if(readId.charAt(i) == READ_ID_DELIM)
                ++delimCount;

            if(delimCount == 3)
                return i + 1;
        }

        return -1;
    }
}
