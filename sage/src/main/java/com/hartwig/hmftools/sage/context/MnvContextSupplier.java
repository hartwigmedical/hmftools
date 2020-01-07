package com.hartwig.hmftools.sage.context;

import java.io.File;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.phase.MnvCandidate;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.sam.SamSlicer;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.select.PositionSelector;
import com.hartwig.hmftools.sage.select.TierSelector;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class MnvContextSupplier {

    private static final Logger LOGGER = LogManager.getLogger(MnvContextSupplier.class);

    private final String sample;
    private final String bamFile;
    private final SageConfig config;
    private final SamSlicerFactory samSlicerFactory;

    private final TierSelector tierSelector;
    private final IndexedFastaSequenceFile refGenome;

    public MnvContextSupplier(@NotNull final SageConfig config, @NotNull final String sample, @NotNull final String bamFile,
            @NotNull final SamSlicerFactory samSlicerFactory, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions, final IndexedFastaSequenceFile refGenome) {
        this.config = config;
        this.sample = sample;
        this.bamFile = bamFile;
        this.samSlicerFactory = samSlicerFactory;
        this.tierSelector = new TierSelector(panelRegions, hotspots);
        this.refGenome = refGenome;
    }

    @NotNull
    public List<AltContext> get(@NotNull final MnvCandidate mnv) {

        final TumorRefContextCandidates candidates = new TumorRefContextCandidates(sample);
        final List<AltContext> altContexts = Lists.newArrayList();
        final PositionSelector<AltContext> consumerSelector = new PositionSelector<>(altContexts);

        GenomeRegion bounds = GenomeRegions.create(mnv.chromosome(), mnv.position(), mnv.end());
        RefSequence refSequence = new RefSequence(mnv, refGenome);
        RefContextConsumer refContextConsumer = new RefContextConsumer(true, config, bounds, refSequence, candidates);
        final SamSlicer slicer = samSlicerFactory.create(bounds);

        try (final SamReader tumorReader = SamReaderFactory.makeDefault().open(new File(bamFile))) {

            slicer.slice(tumorReader, x -> refContextConsumer.processTargeted(mnv, x));

            // Add all valid alt contexts
            candidates.refContexts().stream().flatMap(x -> x.alts().stream()).filter(this::altSupportPredicate).forEach(x -> {
                x.setPrimaryReadCounterFromInterim();
                altContexts.add(x);
            });

            slicer.slice(tumorReader, samRecord -> {
                final IndexedBases refBases = refSequence.alignment(samRecord);
                consumerSelector.select(samRecord.getAlignmentStart(),
                        samRecord.getAlignmentEnd(),
                        x -> x.primaryReadContext().accept(x.rawDepth() < config.maxReadDepth(), samRecord, config, refBases));

            });

        } catch (Exception e) {
            throw new CompletionException(e);
        }

        return altContexts.stream().peek(x -> x.localPhaseSet(mnv.lps())).filter(this::qualPredicate).collect(Collectors.toList());

    }

    private boolean altSupportPredicate(@NotNull final AltContext altContext) {
        return altContext.rawAltSupport() >= config.filter().hardMinTumorAltSupport() || tierSelector.isHotspot(altContext);
    }

    private boolean qualPredicate(@NotNull final AltContext altContext) {
        return altContext.primaryReadContext().tumorQuality() >= config.filter().hardMinTumorQual() || tierSelector.isHotspot(altContext);
    }

}
