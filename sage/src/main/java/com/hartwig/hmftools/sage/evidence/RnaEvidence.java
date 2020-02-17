package com.hartwig.hmftools.sage.evidence;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefContextCandidates;
import com.hartwig.hmftools.sage.context.RefContextConsumer;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RnaEvidence {

    private static final Logger LOGGER = LogManager.getLogger(RnaEvidence.class);

    private final SageConfig config;
    private final SamSlicerFactory samSlicerFactory;

    public RnaEvidence(@NotNull final SageConfig config, final SamSlicerFactory samSlicerFactory) {
        this.config = config;
        this.samSlicerFactory = samSlicerFactory;
    }

    @NotNull
    public List<RefContext> get(@NotNull final RefSequence refSequence, @NotNull final GenomeRegion bounds,
            @NotNull final RefContextCandidates candidates) {
        if (config.rnaEnabled()) {
            final RefContextConsumer refContextConsumer = new RefContextConsumer(false, config, bounds, refSequence, candidates);
            return get(bounds, refContextConsumer, candidates);
        }

        return Collections.emptyList();
    }

    @NotNull
    private List<RefContext> get(@NotNull final GenomeRegion bounds, @NotNull final Consumer<SAMRecord> recordConsumer,
            @NotNull final RefContextCandidates candidates) {

        final SamSlicer slicer = samSlicerFactory.create(bounds);
        try (final SamReader bamReader = SamReaderFactory.makeDefault().open(new File(config.rnaBam()))) {
            slicer.slice(bamReader, recordConsumer);
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

}
