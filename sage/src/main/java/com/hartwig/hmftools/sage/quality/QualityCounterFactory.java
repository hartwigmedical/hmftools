package com.hartwig.hmftools.sage.quality;

import java.io.File;
import java.util.Collection;
import java.util.concurrent.CompletionException;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.SageApplication;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

class QualityCounterFactory {

    private static final Logger LOGGER = LogManager.getLogger(SageApplication.class);

    private final String bamFile;
    private final ReferenceSequenceFile refGenome;
    private final int maxAltCount;

    public QualityCounterFactory(final String bamFile, final ReferenceSequenceFile refGenome, final int maxAltCount) {
        this.bamFile = bamFile;
        this.refGenome = refGenome;
        this.maxAltCount = maxAltCount;
    }

    @NotNull
    public Collection<QualityCounter> regionCount(@NotNull final GenomeRegion bounds) {
        LOGGER.debug("Processing bqr region {}", bounds);

        final RefSequence refSequence = new RefSequence(bounds, refGenome);
        final QualityCounterCigarHandler counter = new QualityCounterCigarHandler(refSequence, bounds, maxAltCount);
        final SamSlicer slicer = new SamSlicer(10, bounds);
        try (final SamReader tumorReader = SamReaderFactory.makeDefault()
                .referenceSource(new ReferenceSource(refGenome))
                .open(new File(bamFile))) {

            // First parse
            slicer.slice(tumorReader, counter::processRecord);

        } catch (Exception e) {
            throw new CompletionException(e);
        }

        return counter.counts();
    }

}
