package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.select.PositionSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class MnvRefContextSupplier {

    private static final Logger LOGGER = LogManager.getLogger(MnvRefContextSupplier.class);

    private final SageConfig config;
    private final IndexedFastaSequenceFile refGenome;
    private final int minQuality;
    private final SageConfig sageConfig;

    public MnvRefContextSupplier(@NotNull final SageConfig config, @NotNull final IndexedFastaSequenceFile refGenome) {
        this.minQuality = config.minMapQuality();
        this.sageConfig = config;
        this.refGenome = refGenome;
        this.config = config;

    }

    @NotNull
    public List<RefContext> get(@NotNull final String sample, @NotNull final String bamFile, @NotNull final VariantHotspot target,
            @NotNull final ReadContext readContext) {

        final NormalRefContextCandidates candidates = new NormalRefContextCandidates(sample);
        RefContext refContext = candidates.add(target.chromosome(), target.position());
        refContext.altContext(target.ref(), target.alt()).setPrimaryReadContext(new ReadContextCounter(target, readContext));

        final GenomeRegion bounds = GenomeRegions.create(target.chromosome(), target.position(), target.end());
        final RefSequence refSequence = new RefSequence(target, refGenome);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(false, config, bounds, refSequence, candidates);
        final SamSlicer slicer = new SamSlicer(config.minMapQuality(), bounds);

        final PositionSelector<AltContext> consumerSelector =
                new PositionSelector<>(candidates.refContexts().stream().flatMap(x -> x.alts().stream()).collect(Collectors.toList()));

        try (final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile))) {
            slicer.slice(tumorReader, samRecord -> {

                refContextConsumer.processTargeted(target, samRecord);
                final IndexedBases refBases = refSequence.alignment(samRecord);

                if (samRecord.getMappingQuality() >= minQuality) {
                    consumerSelector.select(samRecord.getAlignmentStart(),
                            samRecord.getAlignmentEnd(),
                            x -> x.primaryReadContext().accept(x.rawDepth() < sageConfig.maxReadDepth(), samRecord, sageConfig, refBases));
                }

            });
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

}
