package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;

import htsjdk.io.HtsPath;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SamSlicerFactory
{
    private final Map<String,SamSlicerInterface> mSamSlicers;

    private final Map<String,SamReader> mBamReaders;

    public SamSlicerFactory()
    {
        mSamSlicers = Maps.newHashMap();
        mBamReaders = Maps.newHashMap();
    }

    public SamSlicerInterface getSamSlicer(final String sampleId, final List<ChrBaseRegion> regions, boolean keepSupplementaries)
    {
        if(!mBamReaders.isEmpty())
        {
            SamReader bamReader = mBamReaders.get(sampleId);
            return new SamSlicer(bamReader, 0, regions, keepSupplementaries);
        }

        return mSamSlicers.get(sampleId);
    }

    public void addSamSlicer(final String sampleId, final SamSlicerInterface samSlicer)
    {
        mSamSlicers.put(sampleId, samSlicer);
    }

    public void buildBamReaders(
            final List<String> tumorIds, final List<String> tumorBams, final SageConfig config, final IndexedFastaSequenceFile refGenome)
    {
        List<String> allSamples = Lists.newArrayList(tumorIds);
        allSamples.addAll(config.ReferenceIds);

        List<String> allBams = Lists.newArrayList(tumorBams);
        allBams.addAll(config.ReferenceBams);

        for(int i = 0; i < allSamples.size(); i++)
        {
            final String sample = allSamples.get(i);
            final String bamFile = allBams.get(i);

            HtsPath path = new HtsPath(bamFile);

            SamReader bamReader = SamReaderFactory.makeDefault()
                    .validationStringency(config.BamStringency)
                    .referenceSource(new ReferenceSource(refGenome))
                    .open(SamInputResource.of(path.getURI()));

            mBamReaders.put(sample, bamReader);
        }
    }

    public void close()
    {
        try
        {
            for(SamReader bamReader : mBamReaders.values())
            {
                bamReader.close();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to close bam file: {}", e.toString());
        }
    }
}
