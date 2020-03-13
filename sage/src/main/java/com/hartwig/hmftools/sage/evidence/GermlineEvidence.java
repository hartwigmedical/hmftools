package com.hartwig.hmftools.sage.evidence;

import java.io.File;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContextConsumer;
import com.hartwig.hmftools.sage.context.RefContextFactory;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.HotspotSelector;
import com.hartwig.hmftools.sage.select.SamRecordSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class GermlineEvidence {

    private static final Logger LOGGER = LogManager.getLogger(GermlineEvidence.class);

    private final SageConfig config;
    private final List<VariantHotspot> hotspots;
    private final ReferenceSequenceFile refGenome;
    private final SamSlicerFactory samSlicerFactory;

    public GermlineEvidence(@NotNull final SageConfig config, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final SamSlicerFactory samSlicerFactory, @NotNull final ReferenceSequenceFile refGenome) {
        this.config = config;
        this.samSlicerFactory = samSlicerFactory;
        this.hotspots = hotspots;
        this.refGenome = refGenome;
    }

    @NotNull
    public List<AltContext> get(@NotNull final String sample, @NotNull final String bamFile, @NotNull final RefSequence refSequence,
            @NotNull final GenomeRegion bounds) {

        if (bounds.start() == 1) {
            LOGGER.info("Beginning processing of {} chromosome {} ", sample, bounds.chromosome());
        }

        LOGGER.debug("Variant candidates {} position {}:{}", sample, bounds.chromosome(), bounds.start());

        final HotspotSelector hotspotSelector = new HotspotSelector(hotspots);
        final RefContextFactory candidates = new RefContextFactory(config, hotspotSelector, sample);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(config, bounds, refSequence, candidates);
        return get(bamFile, bounds, refContextConsumer, candidates, hotspotSelector);
    }

    @NotNull
    private List<AltContext> get(@NotNull final String bamFile, @NotNull final GenomeRegion bounds,
            @NotNull final Consumer<SAMRecord> recordConsumer, @NotNull final RefContextFactory candidates, @NotNull final HotspotSelector hotspotSelector) {
        final List<AltContext> altContexts = Lists.newArrayList();
        final SamRecordSelector<AltContext> consumerSelector = new SamRecordSelector<>(config.maxSkippedReferenceRegions(), altContexts);


        final SamSlicer slicer = samSlicerFactory.create(bounds);
        try (final SamReader tumorReader = SamReaderFactory.makeDefault()
                .referenceSource(new ReferenceSource(refGenome))
                .open(new File(bamFile))) {

            // First parse
            slicer.slice(tumorReader, recordConsumer);

            // Add all valid alt contexts
            altContexts.addAll(candidates.altContexts());

            // Second parse
            slicer.slice(tumorReader, samRecord -> {
                consumerSelector.select(samRecord,
                        x -> x.primaryReadContext().accept(samRecord, config));
            });

        } catch (Exception e) {
            throw new CompletionException(e);
        }

        return altContexts.stream().filter(x -> qualPredicate(hotspotSelector, x)).collect(Collectors.toList());
    }

    private boolean qualPredicate(@NotNull final HotspotSelector tierSelector, @NotNull final AltContext altContext) {
        return altContext.primaryReadContext().tumorQuality() >= config.filter().hardMinTumorQual() || tierSelector.isHotspot(altContext);
    }

}
