package com.hartwig.hmftools.sage.task;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CompletionException;
import java.util.function.Supplier;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.count.BaseDetails;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class CandidateSupplier implements Supplier<List<BaseDetails>> {

    private static final Logger LOGGER = LogManager.getLogger(CandidateSupplier.class);

    private final String sample;
    private final GenomeRegion region;
    private final String bamFile;
    private final CandidateConsumer samConsumer;

    public CandidateSupplier(final String sample, final GenomeRegion region, final String bamFileLocation,
            final IndexedFastaSequenceFile refGenome) {
        this.sample = sample;
        this.region = region;
        this.bamFile = bamFileLocation;
        samConsumer = new CandidateConsumer(13, region, refGenome);
    }

    public CandidateSupplier(final String sample, final GenomeRegion region, final String bamFileLocation,
            final IndexedFastaSequenceFile refGenome, Set<Long> hotspots) {
        this.sample = sample;
        this.region = region;
        this.bamFile = bamFileLocation;
        samConsumer = new CandidateConsumer(13, region, refGenome, hotspots);
    }

    @Override
    public List<BaseDetails> get() {
        try {
            return candidates();
        } catch (IOException e) {
            throw new CompletionException(e);
        }
    }

    @NotNull
    private List<BaseDetails> candidates() throws IOException {

        LOGGER.info("Candidates " + sample + "  from " + region.start());

        SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile));
        SAMSlicer slicer = new SAMSlicer(0, Lists.newArrayList(region));
        slicer.slice(tumorReader, samConsumer);
        tumorReader.close();

        return samConsumer.bases();
    }

}
