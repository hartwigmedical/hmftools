package com.hartwig.hmftools.isofox.refdata;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;

import java.io.BufferedWriter;

public class RefDataWriter
{
    private final RefDataConfig mConfig;
    private BufferedWriter mExpRateWriter;
    private BufferedWriter mGcRatioWriter;

    public RefDataWriter(final RefDataConfig config)
    {
        mConfig = config;

        if(mConfig.GenerateExpectedCounts)
        {
            mExpRateWriter = ExpectedRatesGenerator.createWriter(mConfig);
        }

        if(mConfig.GenerateGcRatios)
        {
            mGcRatioWriter = null; // ExpectedRatesGenerator.createWriter(mConfig);
        }

    }

    public BufferedWriter getExpRatesWriter() { return mExpRateWriter;}
    public BufferedWriter getReadGcRatioWriter() { return mGcRatioWriter; }

    public void close()
    {
        closeBufferedWriter(mExpRateWriter);
        closeBufferedWriter(mGcRatioWriter);
    }
}
