package com.hartwig.hmftools.sage.phase;

import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.NormalRefContextCandidates;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefSequence;
import com.hartwig.hmftools.sage.evidence.NormalEvidence;
import com.hartwig.hmftools.sage.evidence.PrimaryEvidence;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.sam.SamSlicerFactory;
import com.hartwig.hmftools.sage.sam.SamSlicerFactoryChromImpl;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class MnvFactory {

    private static final Logger LOGGER = LogManager.getLogger(MnvFactory.class);

    private final SageConfig config;
    private final SageVariantFactory sageVariantFactory;
    private final IndexedFastaSequenceFile refGenome;
    private final PrimaryEvidence primaryEvidence;
    private final NormalEvidence normalEvidence;

    MnvFactory(@NotNull final SageConfig config, @NotNull final SageVariantFactory sageVariantFactory,
            @NotNull final IndexedFastaSequenceFile refGenome, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions) {
        final SamSlicerFactory slicerFactory = new SamSlicerFactoryChromImpl(config, panelRegions);
        this.config = config;
        this.sageVariantFactory = sageVariantFactory;
        this.refGenome = refGenome;
        this.primaryEvidence = new PrimaryEvidence(config, hotspots, panelRegions, slicerFactory);
        this.normalEvidence = new NormalEvidence(config, slicerFactory);
    }

    @Nullable
    public SageVariant mnv(int lps, @NotNull final VariantHotspot mnv) {

        final RefSequence refSequence = new RefSequence(mnv, refGenome);

        final List<AltContext> tumorAltContexts = Lists.newArrayList();
        for (int sampleNumber = 0; sampleNumber < config.tumor().size(); sampleNumber++) {

            String sample = config.tumor().get(sampleNumber);
            String bamFile = config.tumorBam().get(sampleNumber);

            final List<AltContext> sampleMnv = primaryEvidence.get(sample, bamFile, refSequence, mnv);
            if (sampleMnv.isEmpty()) {
                return null;
            }

            tumorAltContexts.add(sampleMnv.get(0));
        }

        final ReadContext primaryReadContext = tumorAltContexts.stream()
                .map(AltContext::primaryReadContext)
                .sorted(Comparator.comparingInt(ReadContextCounter::altSupport).reversed())
                .map(ReadContextCounter::readContext)
                .findFirst()
                .orElse(tumorAltContexts.get(0).primaryReadContext().readContext());

        final NormalRefContextCandidates candidates = new NormalRefContextCandidates(config.reference());
        RefContext refContext = candidates.add(mnv.chromosome(), mnv.position());
        refContext.altContext(mnv.ref(), mnv.alt()).setPrimaryReadContext(new ReadContextCounter(mnv, primaryReadContext));

        final List<RefContext> normalRefContexts = normalEvidence.get(config.referenceBam(), refSequence, mnv, candidates);
        final AltContext normalAltContext =
                normalRefContexts.stream().flatMap(x -> x.alts().stream()).findFirst().orElse(new AltContext(config.reference(), mnv));

        final SageVariant result = sageVariantFactory.create(normalAltContext, tumorAltContexts);
        result.localPhaseSet(lps);
        result.synthetic(true);
        return result;
    }

    @NotNull
    public VariantHotspot merge(@NotNull final AltContext left, @NotNull final AltContext right) {
        int mnvLength = (int) (right.position() - left.position() + 1);
        int additionalLength = mnvLength - left.alt().length();

        try {
            final String alt = left.alt() + right.primaryReadContext().readContext().mnvAdditionalAlt(additionalLength);
            final String ref = refGenome.getSubsequenceAt(left.chromosome(), left.position(), right.position()).getBaseString();
            return ImmutableVariantHotspotImpl.builder().from(left).ref(ref).alt(alt).build();
        } catch (Exception e) {
            LOGGER.error("Unable to merge {}:{} with {}", left.chromosome(), left.position(), right.position());
            throw e;
        }
    }

}
