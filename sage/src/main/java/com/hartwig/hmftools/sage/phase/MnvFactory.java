package com.hartwig.hmftools.sage.phase;

import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.MnvAltContextSupplier;
import com.hartwig.hmftools.sage.context.MnvRefContextSupplier;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class MnvFactory {

    private static final Logger LOGGER = LogManager.getLogger(MnvFactory.class);

    private final SageConfig config;
    private final SageVariantFactory sageVariantFactory;
    private final IndexedFastaSequenceFile refGenome;
    private final MnvAltContextSupplier altContextSupplier;
    private final MnvRefContextSupplier refContextSupplier;

    public MnvFactory(@NotNull final SageConfig config, @NotNull final SageVariantFactory sageVariantFactory,
            @NotNull final IndexedFastaSequenceFile refGenome, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions) {
        this.config = config;
        this.sageVariantFactory = sageVariantFactory;
        this.refGenome = refGenome;
        this.altContextSupplier = new MnvAltContextSupplier(config, hotspots, panelRegions, refGenome);
        this.refContextSupplier = new MnvRefContextSupplier(config, refGenome);
    }

    @NotNull
    public SageVariant mnv(int lps, @NotNull final VariantHotspot mnv) {

        final List<AltContext> tumorAltContexts = Lists.newArrayList();
        for (int sampleNumber = 0; sampleNumber < config.tumor().size(); sampleNumber++) {

            String sample = config.tumor().get(sampleNumber);
            String bamFile = config.tumorBam().get(sampleNumber);

            final List<AltContext> sampleMnv = altContextSupplier.get(sample, bamFile, mnv);
            AltContext altContext = sampleMnv.isEmpty() ? new AltContext(sample, mnv) : sampleMnv.get(0);
            tumorAltContexts.add(altContext);
        }

        final ReadContext primaryReadContext = tumorAltContexts.stream()
                .map(AltContext::primaryReadContext)
                .sorted(Comparator.comparingInt(ReadContextCounter::altSupport).reversed())
                .map(ReadContextCounter::readContext)
                .findFirst()
                .orElse(tumorAltContexts.get(0).primaryReadContext().readContext());

        final List<RefContext> normalRefContexts =
                refContextSupplier.get(config.reference(), config.referenceBam(), mnv, primaryReadContext);
        final AltContext normalAltContext =
                normalRefContexts.stream().flatMap(x -> x.alts().stream()).findFirst().orElse(new AltContext(config.reference(), mnv));

        final SageVariant result = sageVariantFactory.create(normalAltContext, tumorAltContexts);
        result.localPhaseSet(lps);
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
