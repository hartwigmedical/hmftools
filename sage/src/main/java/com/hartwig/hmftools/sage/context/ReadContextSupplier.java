package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ReadContextSupplier implements Supplier<List<AltContext>>, Consumer<SAMRecord> {

    private static final Logger LOGGER = LogManager.getLogger(ReadContextSupplier.class);

    private final String sample;
    private final GenomeRegion bounds;
    private final String bamFile;
    private final GenomePositionSelector<ReadContextCounter> consumerSelector;
    private final List<AltContext> altContexts;
    private final int minQuality;

    public ReadContextSupplier(final int minQuality, final String sample, @NotNull final GenomeRegion bounds, @NotNull final String bamFile,
            @NotNull final List<AltContext> altContexts) {

        this.minQuality = minQuality;
        this.sample = sample;
        this.bounds = bounds;
        this.bamFile = bamFile;
        this.altContexts = altContexts;
        consumerSelector =
                GenomePositionSelectorFactory.create(altContexts.stream().map(AltContext::primaryReadContext).collect(Collectors.toList()));

    }

    @Override
    public List<AltContext> get() {

        LOGGER.info("Read Contexts " + sample + "  from " + bounds.start());
        altContexts.forEach(x -> x.primaryReadContext().reset());

        try {
            SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile));
            SAMSlicer slicer = new SAMSlicer(minQuality, Lists.newArrayList(bounds));
            slicer.slice(tumorReader, this);
            tumorReader.close();
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return altContexts;
    }

    @Override
    public void accept(final SAMRecord samRecord) {
        final GenomeRegion samRegion =
                GenomeRegions.create(samRecord.getContig(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());
        consumerSelector.select(samRegion, x -> x.accept(samRecord));
    }
}
