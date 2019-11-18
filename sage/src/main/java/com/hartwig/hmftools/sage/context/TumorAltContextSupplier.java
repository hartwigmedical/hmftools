package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.PositionSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorAltContextSupplier implements Supplier<List<AltContext>> {

    private static final Logger LOGGER = LogManager.getLogger(TumorAltContextSupplier.class);

    private final String sample;
    private final SageConfig config;
    private final String bamFile;
    private final List<AltContext> altContexts = Lists.newArrayList();
    private final PositionSelector<AltContext> consumerSelector;
    private final TumorRefContextCandidates candidates;
    private final RefContextConsumer refContextConsumer;
    private final SamSlicerFactory samSlicerFactory;

    private final GenomeRegion bounds;

    public TumorAltContextSupplier(final SageConfig config, final String sample, @NotNull final GenomeRegion bounds,
            @NotNull final String bamFile, @NotNull final RefSequence refGenome, @NotNull final SamSlicerFactory samSlicerFactory) {
        this.config = config;
        this.sample = sample;
        this.bamFile = bamFile;
        this.samSlicerFactory = samSlicerFactory;
        this.consumerSelector = new PositionSelector<>(altContexts);
        this.candidates = new TumorRefContextCandidates(sample);
        this.bounds = bounds;
        this.refContextConsumer = new RefContextConsumer(true, config, bounds, refGenome, this.candidates);

    }

    private void processFirstPass(final SAMRecord samRecord) {
        refContextConsumer.accept(samRecord);
    }

    private void processSecondPass(final SAMRecord samRecord) {
        consumerSelector.select(samRecord.getAlignmentStart(),
                samRecord.getAlignmentEnd(),
                x -> x.primaryReadContext().accept(x.readDepth() < config.maxReadDepth(), samRecord, config));
    }

    @Override
    public List<AltContext> get() {

        final SamSlicer slicer = samSlicerFactory.create(bounds);

        if (bounds.start() == 1) {
            LOGGER.info("Beginning processing of {} chromosome {} ", sample, bounds.chromosome());
        }

        LOGGER.info("Tumor candidates {} position {}:{}", sample, bounds.chromosome(), bounds.start());

        try (final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile))) {

            slicer.slice(tumorReader, this::processFirstPass);

            // Add all valid alt contexts
            for (final RefContext refContext : candidates.refContexts()) {
                for (final AltContext altContext : refContext.alts()) {
                    if (altContext.altSupport() >= config.filter().hardMinTumorAltSupport()) {
                        altContext.setPrimaryReadCounterFromInterim();
                        altContexts.add(altContext);
                    }
                }
            }

            slicer.slice(tumorReader, this::processSecondPass);

        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return hardQualFilter(altContexts);
    }

    @NotNull
    private List<AltContext> hardQualFilter(@NotNull final List<AltContext> altContexts) {
        return altContexts.stream()
                .filter(x -> x.primaryReadContext().quality() >= config.filter().hardMinTumorQual())
                .collect(Collectors.toList());
    }

}
