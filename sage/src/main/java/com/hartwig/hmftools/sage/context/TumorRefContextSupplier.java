package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorRefContextSupplier implements Supplier<List<RefContext>> {

    private static final Logger LOGGER = LogManager.getLogger(TumorRefContextSupplier.class);

    private final String sample;
    private final String bamFile;
    private final GenomeRegion bounds;
    private final TumorRefContextCandidates candidates;
    private final RefContextConsumer refContextConsumer;

    public TumorRefContextSupplier(final int minQuality, final String sample, @NotNull final GenomeRegion bounds,
            @NotNull final String bamFile, @NotNull final RefSequence refGenome) {
        this.sample = sample;
        this.bounds = bounds;
        this.bamFile = bamFile;
        this.candidates = new TumorRefContextCandidates(sample);
        this.refContextConsumer = new RefContextConsumer(true, minQuality, bounds, refGenome, this.candidates);
    }

    @Override
    public List<RefContext> get() {

        if (bounds.start() == 1) {
            LOGGER.info("Beginning processing of {} chromosome {} ", sample, bounds.chromosome());
        }

        LOGGER.info("Tumor candidates {} position {}:{}", sample, bounds.chromosome(), bounds.start());

        try {
            SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile));
            SageSamSlicer slicer = new SageSamSlicer(0, Lists.newArrayList(bounds));
            slicer.slice(tumorReader, refContextConsumer);
            tumorReader.close();
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

    @NotNull
    public List<RefContext> refContexts() {
        return candidates.refContexts();
    }

}
