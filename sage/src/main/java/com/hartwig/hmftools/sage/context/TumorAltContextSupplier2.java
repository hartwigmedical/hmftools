package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.SageConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

@Deprecated
public class TumorAltContextSupplier2 implements Supplier<List<AltContext>>, Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(TumorAltContextSupplier2.class);

    private final String sample;
    private final SageConfig config;
    private final String bamFile;
    private final ArrayDeque<SAMRecord> samRecordQueue;
    private final List<AltContext> altContexts = Lists.newArrayList();
    private final ContextSelector<AltContext> consumerSelector;
    private final TumorRefContextCandidates candidates;
    private final RefContextConsumer refContextConsumer;

    private final GenomeRegion bounds;

    public TumorAltContextSupplier2(final SageConfig config, final String sample, @NotNull final GenomeRegion bounds,
            @NotNull final String bamFile, @NotNull final RefSequence refGenome) {
        this.config = config;
        this.sample = sample;
        this.bamFile = bamFile;
        this.consumerSelector = new ContextSelector<>(altContexts);
        this.samRecordQueue = new ArrayDeque<>(1_000_000);
        this.candidates = new TumorRefContextCandidates(sample);
        this.bounds = bounds;
        this.refContextConsumer = new RefContextConsumer(true, config, bounds, refGenome, this.candidates);

    }

    @Override
    public void accept(final SAMRecord samRecord) {
        samRecordQueue.addLast(samRecord);
        refContextConsumer.accept(samRecord);
    }

    @Override
    public List<AltContext> get() {

        if (bounds.start() == 1) {
            LOGGER.info("Beginning processing of {} chromosome {} ", sample, bounds.chromosome());
        }

        LOGGER.info("Tumor candidates {} position {}:{}", sample, bounds.chromosome(), bounds.start());

        try {
            SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile));
            SageSamSlicer slicer = new SageSamSlicer(0, Lists.newArrayList(bounds));
            slicer.slice(tumorReader, this);
            tumorReader.close();
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        // Add all valid alt contexts to
        for (final RefContext refContext : candidates.refContexts()) {
            for (final AltContext altContext : refContext.alts()) {
                if (altContext.altSupport() >= config.minTumorAltSupport()) {
                    altContext.setPrimaryReadCounterFromInterim();
                    altContexts.add(altContext);
                }
            }
        }

        // Replay all sam records a second time for the accurate read context counts
        while (!samRecordQueue.isEmpty()) {
            SAMRecord record = samRecordQueue.removeFirst();
            consumerSelector.select(record.getAlignmentStart(), record.getAlignmentEnd(), x -> x.primaryReadContext().accept(record));
        }

        return altContexts;
    }

}
