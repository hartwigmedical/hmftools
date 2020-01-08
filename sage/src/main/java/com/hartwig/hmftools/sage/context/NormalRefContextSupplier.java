package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.PositionSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class NormalRefContextSupplier implements Supplier<List<RefContext>>, Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(NormalRefContextSupplier.class);

    private final GenomeRegion bounds;
    private final RefContextCandidates candidates;
    private final String bamFile;
    private final RefContextConsumer refContextConsumer;
    private final PositionSelector<AltContext> consumerSelector;
    private final int minQuality;
    private final SageConfig sageConfig;
    private final SamSlicerFactory samSlicerFactory;
    private final RefSequence refSequence;

    public NormalRefContextSupplier(final SageConfig config, @NotNull final GenomeRegion bounds, @NotNull final String bamFile,
            @NotNull final RefSequence refSequence, @NotNull final RefContextCandidates candidates,
            @NotNull final SamSlicerFactory samSlicerFactory) {
        this.minQuality = config.minMapQuality();
        this.bounds = bounds;
        this.candidates = candidates;
        this.bamFile = bamFile;
        this.sageConfig = config;
        this.samSlicerFactory = samSlicerFactory;
        this.refSequence = refSequence;
        refContextConsumer = new RefContextConsumer(false, config, bounds, refSequence, candidates);
        consumerSelector =
                new PositionSelector<>(candidates.refContexts().stream().flatMap(x -> x.alts().stream()).collect(Collectors.toList()));
    }

    @Override
    public List<RefContext> get() {

        LOGGER.info("Normal candidates position {}:{}", bounds.chromosome(), bounds.start());

        try (final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile))) {
            samSlicerFactory.create(bounds).slice(tumorReader, this);
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

    @Override
    public void accept(final SAMRecord samRecord) {
        refContextConsumer.accept(samRecord);
        final IndexedBases refBases = refSequence.alignment(samRecord);

        if (samRecord.getMappingQuality() >= minQuality) {
            consumerSelector.select(samRecord.getAlignmentStart(),
                    samRecord.getAlignmentEnd(),
                    x -> x.primaryReadContext().accept(x.rawDepth() < sageConfig.maxReadDepth(), samRecord, sageConfig, refBases));
        }
    }
}
