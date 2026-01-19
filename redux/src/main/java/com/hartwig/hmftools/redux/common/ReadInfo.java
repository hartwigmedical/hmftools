package com.hartwig.hmftools.redux.common;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.redux.duplicate.FragmentCoords;
import com.hartwig.hmftools.redux.duplicate.UmiConfig;

import htsjdk.samtools.SAMRecord;

public class ReadInfo
{
    private final SAMRecord mRead;
    private final FragmentCoords mCoordinates;
    private String mUmi;

    public ReadInfo(final SAMRecord read, final FragmentCoords fragCoords)
    {
        mRead = read;
        mCoordinates = fragCoords;
        mUmi = null;
    }

    public String id() { return mRead.getReadName(); }
    public SAMRecord read() { return mRead; }

    public FragmentCoords fragCoordinates() { return mCoordinates; }

    public String umi() { return mUmi; }

    public String getOrExtractUmi(final UmiConfig umiConfig)
    {
        if(mUmi == null)
            mUmi = umiConfig.extractUmiId(mRead.getReadName());

        return mUmi;
    }

    public String toString()
    {
        return String.format("read(%s:%d %s) coords(%s)",
                mRead.getReferenceName(), mRead.getAlignmentStart(), id(), mCoordinates);
    }

    public static String readToString(final SAMRecord read) { return SamRecordUtils.readToString(read); }
}
