package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.util.Collection;
import java.util.concurrent.CompletionException;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

class QualityCounterFactory
{
    private final String mBamFile;
    private final ReferenceSequenceFile mRefGenome;
    private final SageConfig mConfig;

    public QualityCounterFactory(final SageConfig config, final String bamFile, final ReferenceSequenceFile refGenome)
    {
        mBamFile = bamFile;
        mRefGenome = refGenome;
        mConfig = config;
    }

    @NotNull
    public Collection<QualityCounter> regionCount(@NotNull final GenomeRegion bounds)
    {
        SG_LOGGER.debug("Processing bqr region {}", bounds);

        final RefSequence refSequence = new RefSequence(bounds, mRefGenome);

        final QualityCounterCigarHandler counter = new QualityCounterCigarHandler(
                refSequence, bounds, mConfig.BaseQualityRecalibration.MaxAltCount);

        final SamSlicer slicer = new SamSlicer(mConfig.MinMapQuality, bounds);

        try(final SamReader tumorReader = SamReaderFactory.makeDefault()
                .validationStringency(mConfig.Stringency)
                .referenceSource(new ReferenceSource(mRefGenome))
                .open(new File(mBamFile)))
        {

            // First parse
            slicer.slice(tumorReader, counter::processRecord);

        } catch(Exception e)
        {
            throw new CompletionException(e);
        }

        return counter.counts();
    }
}
