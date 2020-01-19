package com.hartwig.hmftools.sage.pipeline;

import java.util.Comparator;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
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
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SomaticPipeline implements SageVariantPipeline {

    private static final Logger LOGGER = LogManager.getLogger(SomaticPipeline.class);

    private final SageConfig config;
    private final Executor executor;
    private final List<VariantHotspot> hotspots;
    private final List<GenomeRegion> panelRegions;
    private final PrimaryEvidence primaryEvidence;
    private final NormalEvidence normalEvidence;
    private final IndexedFastaSequenceFile refGenome;

    SomaticPipeline(@NotNull final SageConfig config, @NotNull final Executor executor, @NotNull final List<VariantHotspot> hotspots,
            @NotNull final List<GenomeRegion> panelRegions, final IndexedFastaSequenceFile refGenome) {
        this.config = config;
        this.executor = executor;
        this.refGenome = refGenome;
        final SamSlicerFactory samSlicerFactory = new SamSlicerFactory(config, panelRegions);
        this.hotspots = hotspots;
        this.panelRegions = panelRegions;
        this.primaryEvidence = new PrimaryEvidence(config, hotspots, panelRegions, samSlicerFactory);
        this.normalEvidence = new NormalEvidence(config, samSlicerFactory);
    }

    @NotNull
    public CompletableFuture<List<SageVariant>> variants(@NotNull final GenomeRegion region, @NotNull final RefSequence refSequence) {

        final SageVariantFactory variantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions);
        final SomaticPipelineData somaticPipelineData = new SomaticPipelineData(config.reference(), config.tumor().size(), variantFactory);
        List<String> samples = config.tumor();
        List<String> bams = config.tumorBam();

        final List<CompletableFuture<List<AltContext>>> tumorFutures = Lists.newArrayList();
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final String bam = bams.get(i);

            CompletableFuture<List<AltContext>> candidateFuture =
                    CompletableFuture.supplyAsync(() -> primaryEvidence.get(sample, bam, refSequence, region), executor);

            tumorFutures.add(candidateFuture);
        }

        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));

        final CompletableFuture<List<RefContext>> normalFuture = doneTumor.thenApply(aVoid -> {

            for (int i = 0; i < tumorFutures.size(); i++) {
                CompletableFuture<List<AltContext>> future = tumorFutures.get(i);
                somaticPipelineData.addTumor(i, future.join());
            }

            return normalEvidence.get(config.referenceBam(), refSequence, region, somaticPipelineData.normalCandidates());
        });

        return normalFuture.thenApply(aVoid -> {

            somaticPipelineData.addNormal(normalFuture.join());

            return somaticPipelineData.results();
        });
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

        final SageVariantFactory sageVariantFactory = new SageVariantFactory(config.filter(), hotspots, panelRegions);
        final SageVariant result = sageVariantFactory.create(normalAltContext, tumorAltContexts);
        result.localPhaseSet(lps);
        result.synthetic(true);

        if (result.normal().primaryReadContext().altSupport() != 0) {
            result.filters().add(SageVCF.NORMAL_SUPPORT);
        }

        return result;
    }

    @NotNull
    @Override
    public VariantHotspot combined(@NotNull final SageVariant left, @NotNull final SageVariant right) {
        return combined(refGenome, left.primaryTumor(), right.primaryTumor());
    }
}
