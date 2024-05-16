package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.FileCommon.createBamSlicer;
import static com.hartwig.hmftools.esvee.common.SvConstants.BAM_HEADER_SAMPLE_INDEX_TAG;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.esvee.AssemblyConfig;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReader implements AutoCloseable
{
    private final AssemblyConfig mConfig;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;
    private boolean mCurrentIsReferenceSample;

    public BamReader(final AssemblyConfig config)
    {
        mConfig = config;

        mSamReaders = Lists.newArrayList();
        mCurrentIsReferenceSample = false;

        List<String> combinedBamFiles = mConfig.combinedBamFiles();
        List<String> combinedSampleId = mConfig.combinedSampleIds();

        for(int i = 0; i < combinedSampleId.size(); ++i)
        {
            String bamFile = combinedBamFiles.get(i);

            SamReader samReader = SamReaderFactory.makeDefault()
                            .validationStringency(mConfig.BamStringency)
                            .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFile));

            samReader.getFileHeader().setAttribute(BAM_HEADER_SAMPLE_INDEX_TAG, i);

            mSamReaders.add(samReader);
        }

        mBamSlicer = createBamSlicer();
    }

    public void sliceBam(final String chromosome, int positionStart, int positionEnd, final Consumer<SAMRecord> consumer)
    {
        int bamPosStart = max(positionStart, 1);
        int bamPosEnd = min(positionEnd, mConfig.RefGenomeCoords.length(chromosome));

        if(bamPosStart > bamPosEnd)
            return;

        for(int i = 0; i < mSamReaders.size(); ++i)
        {
            SamReader reader = mSamReaders.get(i);

            mCurrentIsReferenceSample = i >= mConfig.TumorIds.size();
            mBamSlicer.slice(reader, new ChrBaseRegion(chromosome, positionStart, positionEnd), consumer);
        }
    }

    public boolean currentIsReferenceSample() { return mCurrentIsReferenceSample; }

    @Override
    public void close()
    {
        try
        {
            for(SamReader reader : mSamReaders)
            {
                reader.close();
            }
        }
        catch(IOException e) {}
    }
}
