package com.hartwig.hmftools.common.amber;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.SAMConsumer;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilder;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorContaminationEvidence implements Callable<TumorContaminationEvidence> {

    private final String contig;
    private final String bamFile;
    private final SamReaderFactory samReaderFactory;
    private final Map<BaseDepth, ModifiableBaseDepth> evidenceMap;
    private final GenomePositionSelector<ModifiableBaseDepth> selector;
    private final BaseDepthFactory bafFactory;
    private final SAMConsumer supplier;

    public TumorContaminationEvidence(int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig,
            final String bamFile, final SamReaderFactory samReaderFactory, final List<BaseDepth> baseDepths) {
        this.bafFactory = new BaseDepthFactory(minBaseQuality);
        this.contig = contig;
        this.bamFile = bamFile;
        this.samReaderFactory = samReaderFactory;
        final GenomeRegionBuilder builder = new GenomeRegionBuilder(contig, typicalReadDepth);
        baseDepths.forEach(x -> builder.addPosition(x.position()));
        final List<GenomeRegion> bafRegions1 = builder.build();
        final List<ModifiableBaseDepth> tumorRecords = baseDepths.stream().map(BaseDepthFactory::create).collect(Collectors.toList());

        this.evidenceMap = tumorRecords.stream().collect(Collectors.toMap(x -> x, x -> x));
        this.selector = GenomePositionSelectorFactory.create(tumorRecords);
        this.supplier = new SAMConsumer(minMappingQuality, bafRegions1);
    }

    @NotNull
    public String contig() {
        return contig;
    }

    @NotNull
    public List<TumorContamination> evidence() {
        final List<TumorContamination> result = Lists.newArrayList();
        for (final Map.Entry<BaseDepth, ModifiableBaseDepth> entry : evidenceMap.entrySet()) {
            final BaseDepth normal = entry.getKey();
            final BaseDepth tumor = entry.getValue();
            if (tumor.altSupport() != 0) {
                result.add(ImmutableTumorContamination.builder().from(normal).normal(normal).tumor(tumor).build());
            }
        }

        return result;
    }

    @Override
    public TumorContaminationEvidence call() throws Exception {

        try (SamReader reader = samReaderFactory.open(new File(bamFile))) {
            supplier.consume(reader, this::record);
        }

        return this;
    }

    private void record(@NotNull final SAMRecord record) {
        selector.select(asRegion(record), bafEvidence -> bafFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record) {
        return GenomeRegionFactory.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }
}
