package com.hartwig.hmftools.sage.evidence;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefContextCandidates;
import com.hartwig.hmftools.sage.context.RefContextConsumer;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.SamRecordSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class NormalEvidence {

    private static final Logger LOGGER = LogManager.getLogger(NormalEvidence.class);

    private final int minQuality;
    private final SageConfig config;
    private final SageConfig sageConfig;
    private final SamSlicerFactory samSlicerFactory;

    public NormalEvidence(@NotNull final SageConfig config, final SamSlicerFactory samSlicerFactory) {
        this.minQuality = config.minMapQuality();
        this.sageConfig = config;
        this.config = config;
        this.samSlicerFactory = samSlicerFactory;
    }

    @NotNull
    public List<RefContext> get(@NotNull final String bamFile, @NotNull final RefSequence refSequence, @NotNull final GenomeRegion bounds,
            @NotNull final RefContextCandidates candidates) {
        final RefContextConsumer refContextConsumer = new RefContextConsumer(false, config, bounds, refSequence, candidates);
        return get(bamFile, refSequence, bounds, refContextConsumer, candidates);
    }

    @NotNull
    private List<RefContext> get(@NotNull final String bamFile, @NotNull final RefSequence refSequence, @NotNull final GenomeRegion bounds,
            @NotNull final Consumer<SAMRecord> recordConsumer, @NotNull final RefContextCandidates candidates) {

        final SamSlicer slicer = samSlicerFactory.create(bounds);

        final SamRecordSelector<AltContext> consumerSelector =
                new SamRecordSelector<>(candidates.refContexts().stream().flatMap(x -> x.alts().stream()).collect(Collectors.toList()));

        try (final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile))) {
            slicer.slice(tumorReader, samRecord -> {

                recordConsumer.accept(samRecord);
                final IndexedBases refBases = refSequence.alignment(samRecord);

                if (samRecord.getMappingQuality() >= minQuality) {
                    consumerSelector.select(samRecord,
                            x -> x.primaryReadContext().accept(x.rawDepth() < sageConfig.maxReadDepth(), samRecord, sageConfig, refBases));
                }

            });
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

}
