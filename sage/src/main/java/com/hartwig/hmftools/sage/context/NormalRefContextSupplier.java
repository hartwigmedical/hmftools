package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.sam.SimpleSamSlicer;
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
    private final PositionSelector<ReadContextCounter> consumerSelector;
    private final int minQuality;
    private final SageConfig sageConfig;

    public NormalRefContextSupplier(final SageConfig config, @NotNull final GenomeRegion bounds, @NotNull final String bamFile,
            @NotNull final RefSequence refGenome, @NotNull final RefContextCandidates candidates) {
        this.minQuality = config.minMapQuality();
        this.bounds = bounds;
        this.candidates = candidates;
        this.bamFile = bamFile;
        this.sageConfig = config;
        refContextConsumer = new RefContextConsumer(false, config, bounds, refGenome, candidates);
        consumerSelector = new PositionSelector<>(candidates.refContexts()
                .stream()
                .flatMap(x -> x.alts().stream())
                .map(AltContext::primaryReadContext)
                .collect(Collectors.toList()));

    }

    @Override
    public List<RefContext> get() {

        LOGGER.info("Normal candidates position {}:{}", bounds.chromosome(), bounds.start());

        try (final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile))) {
            new SimpleSamSlicer(0, Lists.newArrayList(bounds)).slice(tumorReader, this);
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

    @NotNull
    public List<RefContext> refContexts() {
        return candidates.refContexts();
    }

    @Override
    public void accept(final SAMRecord samRecord) {
        refContextConsumer.accept(samRecord);

        if (samRecord.getMappingQuality() >= minQuality) {
            consumerSelector.select(samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), x -> x.accept(samRecord, sageConfig));
        }
    }
}
