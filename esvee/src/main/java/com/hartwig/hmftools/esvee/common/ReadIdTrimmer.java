package com.hartwig.hmftools.esvee.common;

public class ReadIdTrimmer
{
    private boolean mEnabled;
    private int mReadIdTrimIndex;
    private String mTrimmedBases;

    public ReadIdTrimmer(boolean enabled)
    {
        mEnabled = enabled;
        mReadIdTrimIndex = -1;
        mTrimmedBases = "";
    }

    public boolean enabled() { return mEnabled; }

    public String trim(final String readId)
    {
        if(!mEnabled)
            return readId;

        if(mReadIdTrimIndex < 0)
        {
            mReadIdTrimIndex = findReadIdTrimIndex(readId);

            if(mReadIdTrimIndex < 0)
            {
                // disable if the first read doesn't match the established pattern
                mEnabled = false;
                return readId;
            }

            mTrimmedBases = readId.substring(0, mReadIdTrimIndex);
        }

        if(mReadIdTrimIndex >= readId.length() || readId.charAt(mReadIdTrimIndex - 1) != READ_ID_DELIM)
        {
            // disable if any reads don't match the established pattern
            mEnabled = false;
            return readId;
        }

        return readId.substring(mReadIdTrimIndex);
    }

    public String restore(final String trimmedReadId) { return mEnabled ? mTrimmedBases + trimmedReadId : trimmedReadId; }

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
