package com.hartwig.hmftools.common.amber;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BaseDepthEvidence implements Callable<BaseDepthEvidence> {

    private final String contig;
    private final String bamFile;
    private final SamReaderFactory samReaderFactory;
    private final List<ModifiableBaseDepth> evidence;
    private final GenomePositionSelector<ModifiableBaseDepth> selector;
    private final BaseDepthFactory bafFactory;
    private final SAMSlicer supplier;

    public BaseDepthEvidence(int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig, final String bamFile,
            final SamReaderFactory samReaderFactory, final List<AmberSite> bafRegions) {
        this.bafFactory = new BaseDepthFactory(minBaseQuality);
        this.contig = contig;
        this.bamFile = bamFile;
        this.samReaderFactory = samReaderFactory;
        final GenomeRegionsBuilder builder = new GenomeRegionsBuilder(typicalReadDepth);
        bafRegions.forEach(builder::addPosition);
        final List<GenomeRegion> bafRegions1 = builder.build();

        this.evidence = bafRegions.stream().map(BaseDepthFactory::fromAmberSite).collect(Collectors.toList());
        this.selector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(evidence));
        this.supplier = new SAMSlicer(minMappingQuality, bafRegions1);
    }

    @NotNull
    public String contig() {
        return contig;
    }

    @NotNull
    public List<BaseDepth> evidence() {
        return new ArrayList<>(evidence);
    }

    @Override
    public BaseDepthEvidence call() throws Exception {
        try (SamReader reader = samReaderFactory.open(new File(bamFile))) {
            supplier.slice(reader, this::record);
        }

        return this;
    }

    private void record(@NotNull final SAMRecord record) {
        selector.select(asRegion(record), bafEvidence -> bafFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record) {
        return GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }
}
