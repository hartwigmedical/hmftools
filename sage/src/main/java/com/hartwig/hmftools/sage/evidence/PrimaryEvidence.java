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
import com.hartwig.hmftools.sage.context.RefContextCandidates;
import com.hartwig.hmftools.sage.context.RefContextConsumer;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.context.TumorRefContextCandidates;
import com.hartwig.hmftools.sage.read.IndexedBases;
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

public class PrimaryEvidence {

    private static final Logger LOGGER = LogManager.getLogger(PrimaryEvidence.class);

    private final SageConfig config;
    private final SamSlicerFactory samSlicerFactory;
    private final List<VariantHotspot> hotspots;

    public PrimaryEvidence(@NotNull final SageConfig config, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final SamSlicerFactory samSlicerFactory) {
        this.config = config;
        this.samSlicerFactory = samSlicerFactory;
        this.hotspots = hotspots;
    }

    @NotNull
    public List<AltContext> get(@NotNull final String sample, @NotNull final String bamFile, @NotNull final RefSequence refSequence,
            @NotNull final GenomeRegion bounds) {

        if (bounds.start() == 1) {
            LOGGER.info("Beginning processing of {} chromosome {} ", sample, bounds.chromosome());
        }
        LOGGER.info("Variant candidates {} position {}:{}", sample, bounds.chromosome(), bounds.start());

        final TumorRefContextCandidates candidates = new TumorRefContextCandidates(sample);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(true, config, bounds, refSequence, candidates);
        return get(bamFile, refSequence, bounds, refContextConsumer, candidates);
    }

    @NotNull
    private List<AltContext> get(@NotNull final String bamFile, @NotNull final RefSequence refSequence, @NotNull final GenomeRegion bounds,
            @NotNull final Consumer<SAMRecord> recordConsumer, @NotNull final RefContextCandidates candidates) {
        final List<AltContext> altContexts = Lists.newArrayList();
        final SamRecordSelector<AltContext> consumerSelector = new SamRecordSelector<>(altContexts);
        final HotspotSelector tierSelector = new HotspotSelector(hotspots);

        final SamSlicer slicer = samSlicerFactory.create(bounds);
        try (final SamReader tumorReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(config.refGenome()))
                .open(new File(bamFile))) {

            // First parse
            slicer.slice(tumorReader, recordConsumer);

            // Add all valid alt contexts
            candidates.refContexts()
                    .stream()
                    .flatMap(x -> x.alts().stream())
                    .filter(this::refPredicate)
                    .filter(x -> rawPredicate(tierSelector, x))
                    .forEach(x -> {
                        x.setPrimaryReadCounterFromInterim();
                        altContexts.add(x);
                    });

            // Second parse
            slicer.slice(tumorReader, samRecord -> {
                final IndexedBases refBases = refSequence.alignment(samRecord);
                consumerSelector.select(samRecord,
                        x -> x.primaryReadContext().accept(x.rawDepth() < config.maxReadDepth(), samRecord, config, refBases));

            });

        } catch (Exception e) {
            throw new CompletionException(e);
        }

        return altContexts.stream().filter(x -> qualPredicate(tierSelector, x)).collect(Collectors.toList());
    }

    private boolean refPredicate(@NotNull final AltContext altContext) {
        for (int i = 0; i < altContext.ref().length(); i++) {
            char base = altContext.ref().charAt(i);
            if (base != 'G' && base != 'A' && base != 'T' && base != 'C') {
                return false;
            }
        }

        return true;
    }

    private boolean rawPredicate(@NotNull final HotspotSelector tierSelector, @NotNull final AltContext altContext) {
        return tierSelector.isHotspot(altContext) || altContext.rawSupportAlt() >= config.filter().hardMinTumorRawAltSupport()
                && altContext.rawBaseQualityAlt() >= config.filter().hardMinTumorRawBaseQuality();
    }

    private boolean qualPredicate(@NotNull final HotspotSelector tierSelector, @NotNull final AltContext altContext) {
        return altContext.primaryReadContext().tumorQuality() >= config.filter().hardMinTumorQual() || tierSelector.isHotspot(altContext);
    }

}
