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
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.sam.SimpleSamSlicer;
import com.hartwig.hmftools.sage.select.PositionSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorReadContextSupplier implements Supplier<List<AltContext>>, Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(TumorReadContextSupplier.class);

    private final String sample;
    private final GenomeRegion bounds;
    private final String bamFile;
    private final PositionSelector<ReadContextCounter> consumerSelector;
    private final List<AltContext> altContexts;
    private final int minQuality;

    public TumorReadContextSupplier(final int minQuality, final String sample, @NotNull final GenomeRegion bounds,
            @NotNull final String bamFile, @NotNull final List<AltContext> altContexts) {

        this.minQuality = minQuality;
        this.sample = sample;
        this.bounds = bounds;
        this.bamFile = bamFile;
        this.altContexts = altContexts;
        final List<ReadContextCounter> readContextCounters =
                altContexts.stream().map(AltContext::setPrimaryReadCounterFromInterim).collect(Collectors.toList());
        consumerSelector = new PositionSelector<>(readContextCounters);
    }

    @Override
    public List<AltContext> get() {

        LOGGER.info("Tumor read contexts {} position {}:{}", sample, bounds.chromosome(), bounds.start());
        try {
            SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile));
            SimpleSamSlicer slicer = new SimpleSamSlicer(minQuality, Lists.newArrayList(bounds));
            slicer.slice(tumorReader, this);
            tumorReader.close();
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return altContexts;
    }

    @Override
    public void accept(final SAMRecord samRecord) {
        consumerSelector.select(samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), x -> x.accept(samRecord));
    }
}
