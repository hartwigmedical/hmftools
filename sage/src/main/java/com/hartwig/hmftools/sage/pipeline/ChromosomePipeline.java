package com.hartwig.hmftools.sage.pipeline;

import java.io.File;
import java.io.IOException;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.MitochondrialChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.phase.Phase;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantContextFactory;
import com.hartwig.hmftools.sage.variant.SageVariantTier;
import com.hartwig.hmftools.sage.vcf.SageChromosomeVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomePipeline implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(ChromosomePipeline.class);
    private static final EnumSet<SageVariantTier> PANEL_ONLY_TIERS = EnumSet.of(SageVariantTier.HOTSPOT, SageVariantTier.PANEL);

    private final String chromosome;
    private final SageConfig config;
    private final SageChromosomeVCF sageVCF;
    private final List<CompletableFuture<List<SageVariant>>> regions = Lists.newArrayList();
    private final IndexedFastaSequenceFile refGenome;
    private final SageVariantPipeline sageVariantPipeline;

    public ChromosomePipeline(@NotNull final String chromosome, @NotNull final SageConfig config, @NotNull final Executor executor,
            @NotNull final List<VariantHotspot> hotspots, @NotNull final List<GenomeRegion> panelRegions,
            @NotNull final List<GenomeRegion> highConfidenceRegions, final Map<String, QualityRecalibrationMap> qualityRecalibrationMap)
            throws IOException {
        this.chromosome = chromosome;
        this.config = config;
        this.sageVCF = new SageChromosomeVCF(chromosome, config);
        this.refGenome = new IndexedFastaSequenceFile(new File(config.refGenome()));
        this.sageVariantPipeline =
                new SomaticPipeline(config, executor, refGenome, hotspots, panelRegions, highConfidenceRegions, qualityRecalibrationMap);
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    @NotNull
    public String vcfFilename() {
        return sageVCF.filename();
    }

    public void addAllRegions() {
        int maxPosition = refGenome.getSequence(chromosome).length();
        addAllRegions(maxPosition);
    }

    public void addAllRegions(int maxPosition) {
        // This is for the benefit of MT
        int dynamicSliceSize = maxPosition / Math.min(config.threads(), 4) + 1;

        final int regionSliceSize = Math.min(dynamicSliceSize, config.regionSliceSize());
        for (int i = 0; ; i++) {
            int start = 1 + i * regionSliceSize;
            int end = Math.min(start + regionSliceSize - 1, maxPosition);
            addRegion(start, end);

            if (end >= maxPosition) {
                break;
            }
        }
    }

    public void addRegion(int start, int end) {
        final GenomeRegion region = GenomeRegions.create(chromosome, start, end);
        regions.add(sageVariantPipeline.variants(region));
    }

    @NotNull
    public CompletableFuture<ChromosomePipeline> submit() {

        final Consumer<SageVariant> phasedConsumer = variant -> {
            if (include(variant)) {
                final VariantContext context = SageVariantContextFactory.create(variant);
                sageVCF.write(context);
            }
        };

        final Phase phase = new Phase(config, chromosome, phasedConsumer);

        // Phasing must be done in (positional) order but we can do it eagerly as each new region comes in.
        // It is not necessary to wait for the entire chromosome to be finished to start.
        CompletableFuture<Void> done = CompletableFuture.completedFuture(null);
        Iterator<CompletableFuture<List<SageVariant>>> regionsIterator = regions.iterator();
        while (regionsIterator.hasNext()) {
            CompletableFuture<List<SageVariant>> region = regionsIterator.next();
            done = done.thenCombine(region, (aVoid, sageVariants) -> {

                sageVariants.forEach(phase);
                return null;
            });

            regionsIterator.remove();
        }

        return done.thenApply(aVoid -> {
            phase.flush();
            sageVCF.close();
            LOGGER.info("Processing chromosome {} complete", chromosome);
            return ChromosomePipeline.this;
        });

    }

    private boolean include(@NotNull final SageVariant entry) {

        if (config.panelOnly() && !PANEL_ONLY_TIERS.contains(entry.tier())) {
            return false;
        }

        if (entry.isPassing()) {
            return true;
        }

        if (config.filter().hardFilter()) {
            return false;
        }

        if (entry.tier() == SageVariantTier.HOTSPOT) {
            return true;
        }

        // Its not always 100% transparent whats happening with the mixed germline dedup logic unless we keep all the associated records.
        if (entry.mixedGermlineImpact() > 0) {
            return true;
        }

        if (!entry.isNormalEmpty() && !entry.isTumorEmpty() && !MitochondrialChromosome.contains(entry.chromosome())) {
            final ReadContextCounter normal = entry.normalAltContexts().get(0);
            if (normal.altSupport() > config.filter().hardMaxNormalAltSupport()) {
                return false;
            }
        }

        return true;

    }

    @Override
    public void close() throws IOException {
        refGenome.close();
    }
}
