package com.hartwig.hmftools.esvee.read;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.SvConstants.BAM_HEADER_SAMPLE_ID_TAG;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.esvee.SvConfig;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReader implements AutoCloseable
{
    private final SvConfig mConfig;

    private final List<SamReader> mSamReaders;
    private final BamSlicer mBamSlicer;
    private boolean mCurrentIsReferenceSample;

    public BamReader(final SvConfig config)
    {
        mConfig = config;

        mSamReaders = Lists.newArrayList();
        mCurrentIsReferenceSample = false;

        List<String> combinedBamFiles = mConfig.combinedBamFiles();
        List<String> combinedSampleId = mConfig.combinedSampleIds();

        for(int i = 0; i < combinedSampleId.size(); ++i)
        {
            String sampleId = combinedSampleId.get(i);
            String bamFile = combinedBamFiles.get(i);

            SamReader samReader = SamReaderFactory.makeDefault()
                            .validationStringency(mConfig.BamStringency)
                            .referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(bamFile));

            samReader.getFileHeader().setAttribute(BAM_HEADER_SAMPLE_ID_TAG, sampleId);

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

        for(SamReader reader : mSamReaders)
        {
            String fileSampleId = reader.getFileHeader().getAttribute(BAM_HEADER_SAMPLE_ID_TAG);
            mCurrentIsReferenceSample = mConfig.ReferenceIds.contains(fileSampleId);
            mBamSlicer.slice(reader, new ChrBaseRegion(chromosome, positionStart, positionEnd), consumer);
        }
    }

    public boolean currentIsReferenceSample() { return mCurrentIsReferenceSample; }

    private BamSlicer createBamSlicer()
    {
        // as per SvPrep
        BamSlicer bamSlicer = new BamSlicer(0, false, true, false);
        bamSlicer.setKeepUnmapped();
        bamSlicer.setKeepHardClippedSecondaries();
        return bamSlicer;
    }

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
