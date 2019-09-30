package com.hartwig.hmftools.sage.task;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.Executor;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.count.BaseDetails;
import com.hartwig.hmftools.sage.evidence.ImmutableSampleEvidence;
import com.hartwig.hmftools.sage.evidence.ImmutableVariantEvidence;
import com.hartwig.hmftools.sage.evidence.SampleEvidence;
import com.hartwig.hmftools.sage.evidence.VariantEvidence;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SagePipeline {

    private static final Logger LOGGER = LogManager.getLogger(SagePipeline.class);

    private final GenomeRegion region;
    private final SageConfig config;
    private final Executor executor;
    private final IndexedFastaSequenceFile refGenome;

    public SagePipeline(final GenomeRegion region, final SageConfig config, final Executor executor,
            final IndexedFastaSequenceFile refGenome) {
        this.region = region;
        this.config = config;
        this.executor = executor;
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<VariantEvidence>> submit() {

        LOGGER.info("Starting pipeline of: " + region.start());

        List<String> samples = config.tumor();
        List<String> bams = config.tumorBam();

        final List<CompletableFuture<List<VariantHotspotEvidence>>> tumorFutures = Lists.newArrayList();
        for (int i = 0; i < samples.size(); i++) {
            final String sample = samples.get(i);
            final String bam = bams.get(i);
            tumorFutures.add(CompletableFuture.completedFuture(sample)
                    .thenApplyAsync(unused -> new CandidateSupplier(sample, region, bam, refGenome).get()
                            .stream()
                            .flatMap(x -> x.evidence().stream())
                            .filter(x -> x.altSupport() > 2)
                            .collect(Collectors.toList()), executor));
        }

        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));

        final CompletableFuture<List<BaseDetails>> normalFuture = doneTumor.thenApplyAsync(aVoid -> {
            Set<Long> positions =
                    tumorFutures.stream().flatMap(x -> x.join().stream()).map(VariantHotspotEvidence::position).collect(Collectors.toSet());

            return new CandidateSupplier(config.reference(), region, config.referenceBam(), refGenome, positions).get();
        });

        return normalFuture.thenApplyAsync(aVoid -> {

            LOGGER.info("Creating set " + region.start());

            return combine(normalFuture.join(), tumorFutures.stream().map(CompletableFuture::join).collect(Collectors.toList()));
        });
    }

    @NotNull
    private List<VariantEvidence> combine(@NotNull final List<BaseDetails> normal, @NotNull final List<List<VariantHotspotEvidence>> tumors) {
        final List<VariantEvidence> result = Lists.newArrayList();

        final Set<VariantHotspot> allHotspots = Sets.newHashSet();
        final Map<String, Map<VariantHotspot, VariantHotspotEvidence>> evidencePerSample = Maps.newHashMap();
        for (int i = 0; i < config.tumor().size(); i++) {
            String sample = config.tumor().get(i);

            Map<VariantHotspot, VariantHotspotEvidence> mapEntry = evidenceMap(tumors.get(i));

            evidencePerSample.put(sample, mapEntry);
            allHotspots.addAll(mapEntry.keySet());
        }

        final Map<Long, BaseDetails> normalMap = normal.stream().collect(Collectors.toMap(BaseDetails::position, x -> x));

        for (VariantHotspot hotspot : allHotspots) {
            BaseDetails normalBaseDetails = normalMap.get(hotspot.position());

            final SampleEvidence normalEvidence = createEvidence(config.reference(),
                    hotspot,
                    normalBaseDetails == null ? null : normalBaseDetails.selectOrCreate(hotspot.ref(), hotspot.alt()));
            final List<SampleEvidence> tumorEvidence = Lists.newArrayList();
            for (int i = 0; i < config.tumor().size(); i++) {
                String sample = config.tumor().get(i);
                tumorEvidence.add(createEvidence(sample, hotspot, evidencePerSample.get(sample).get(hotspot)));
            }

            VariantEvidence evidence = ImmutableVariantEvidence.builder()
                    .from(hotspot)
                    .readContext(Strings.EMPTY)
                    .normalEvidence(normalEvidence)
                    .addAllTumorEvidence(tumorEvidence)
                    .build();

            result.add(evidence);
        }

        result.sort(Comparator.comparingLong(GenomePosition::position));

        return result;
    }

    @NotNull
    private Map<VariantHotspot, VariantHotspotEvidence> evidenceMap(@NotNull final List<VariantHotspotEvidence> supplier) {
        return supplier.stream().collect(Collectors.toMap(x -> (VariantHotspot) x, x -> x));
    }

    private SampleEvidence createEvidence(@NotNull final String sample, @NotNull VariantHotspot hotspot,
            @Nullable VariantHotspotEvidence evidence) {

        ImmutableSampleEvidence.Builder builder = ImmutableSampleEvidence.builder().from(hotspot).sample(sample);
        if (evidence == null) {
            builder.readContextFull(0)
                    .readContextPartial(0)
                    .readContextRealigned(0)
                    .readDepth(0)
                    .subprimeReadDepth(0)
                    .refQuality(0)
                    .refSupport(0)
                    .altSupport(0)
                    .altQuality(0);
        } else {
            builder.readContextFull(0)
                    .readContextPartial(0)
                    .readContextRealigned(0)
                    .readDepth(evidence.readDepth())
                    .subprimeReadDepth(evidence.subprimeReadDepth())
                    .refQuality(evidence.refQuality())
                    .refSupport(evidence.refSupport())
                    .altSupport(evidence.altSupport())
                    .altQuality(evidence.altQuality());
        }

        return builder.build();

    }

    //    @NotNull
    //    public CompletableFuture<List<VariantEvidence>> submit() {
    //
    //        LOGGER.info("Starting pipeline of: " + region.start());
    //
    //        List<String> samples = config.tumor();
    //        List<String> bams = config.tumorBam();
    //
    //
    //        final List<CompletableFuture<List<BaseDetails>>> tumorFutures = Lists.newArrayList();
    //        for (int i = 0; i < samples.size(); i++) {
    //            final String sample = samples.get(i);
    //            final String bam = bams.get(i);
    //            tumorFutures.add(CompletableFuture.completedFuture(sample)
    //                    .thenApplyAsync(unused -> new CandidateSupplier(sample, region, bam, refGenome).get(), executor));
    //        }
    //
    //        final CompletableFuture<Void> doneTumor = CompletableFuture.allOf(tumorFutures.toArray(new CompletableFuture[tumorFutures.size()]));
    //
    //        final CompletableFuture<List<BaseDetails>> normalFuture = doneTumor.thenApplyAsync(aVoid -> {
    //            Set<Long> positions =
    //                    tumorFutures.stream().flatMap(x -> x.join().stream()).map(BaseDetails::position).collect(Collectors.toSet());
    //
    //            return new CandidateSupplier(config.reference(), region, config.referenceBam(), refGenome, positions).get();
    //        });
    //
    //        return normalFuture.thenApplyAsync(aVoid -> {
    //
    //            LOGGER.info("Creating set " + region.start());
    //
    //            return combine(normalFuture.join(), tumorFutures.stream().map(CompletableFuture::join).collect(Collectors.toList()));
    //        });
    //    }
}
