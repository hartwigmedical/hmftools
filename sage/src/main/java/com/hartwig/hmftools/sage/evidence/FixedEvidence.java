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
import com.hartwig.hmftools.sage.context.RefContextConsumer;
import com.hartwig.hmftools.sage.context.RefContextFixedFactory;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.SamRecordSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class FixedEvidence {

    private static final Logger LOGGER = LogManager.getLogger(FixedEvidence.class);

    private final int minQuality;
    private final SageConfig sageConfig;
    private final SamSlicerFactory samSlicerFactory;
    private final ReferenceSequenceFile refGenome;

    public FixedEvidence(@NotNull final SageConfig config, @NotNull final SamSlicerFactory samSlicerFactory,
            @NotNull final ReferenceSequenceFile refGenome) {
        this.minQuality = config.minMapQuality();
        this.sageConfig = config;
        this.samSlicerFactory = samSlicerFactory;
        this.refGenome = refGenome;
    }

    @NotNull
    public List<RefContext> get(@NotNull final RefSequence refSequence, @NotNull final GenomeRegion bounds,
            @NotNull final RefContextFixedFactory candidates, @NotNull final String bam) {
        final RefContextConsumer refContextConsumer = new RefContextConsumer(false, sageConfig, bounds, refSequence, candidates);
        return get(bounds, refContextConsumer, candidates, bam);
    }

    @NotNull
    private List<RefContext> get(@NotNull final GenomeRegion bounds, @NotNull final Consumer<SAMRecord> recordConsumer,
            @NotNull final RefContextFixedFactory candidates, @NotNull final String bam) {

        final SamSlicer slicer = samSlicerFactory.create(bounds);

        final SamRecordSelector<AltContext> consumerSelector = new SamRecordSelector<>(sageConfig.maxSkippedReferenceRegions(),
                candidates.refContexts().stream().flatMap(x -> x.alts().stream()).collect(Collectors.toList()));

        try (final SamReader tumorReader = SamReaderFactory.makeDefault()
                .referenceSource(new ReferenceSource(refGenome))
                .open(new File(bam))) {
            slicer.slice(tumorReader, samRecord -> {

                recordConsumer.accept(samRecord);

                if (samRecord.getMappingQuality() >= minQuality) {
                    consumerSelector.select(samRecord,
                            x -> x.primaryReadContext().accept(x.rawDepth() < sageConfig.maxReadDepth(), samRecord, sageConfig));
                }

            });
        } catch (IOException e) {
            throw new CompletionException(e);
        }

        return candidates.refContexts();
    }

}
