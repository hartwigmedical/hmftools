package com.hartwig.hmftools.common.amber;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.hotspot.SAMSupplier;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionBuilder;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TumorBAFEvidence implements Callable<TumorBAFEvidence> {

    private final String contig;
    private final String bamFile;
    private final TumorBAFFactory bafFactory;
    private final List<ModifiableTumorBAF> evidence;
    private final SamReaderFactory samReaderFactory;
    private final GenomePositionSelector<ModifiableTumorBAF> selector;
    private final SAMSupplier supplier;

    public TumorBAFEvidence(int minMappingQuality, int minBaseQuality, final String contig, final String bamFile,
            final SamReaderFactory samReaderFactory, final List<NormalBAF> normalBafs) {
        this.bafFactory = new TumorBAFFactory(minBaseQuality);
        this.contig = contig;
        this.bamFile = bamFile;
        this.samReaderFactory = samReaderFactory;

        final GenomeRegionBuilder builder = new GenomeRegionBuilder(contig, 1000);
        for (NormalBAF bafRegion : normalBafs) {
            builder.addPosition(bafRegion.position());
        }

        final List<GenomeRegion> bafRegions = builder.build();
        this.evidence = normalBafs.stream().map(TumorBAFFactory::create).collect(Collectors.toList());
        this.selector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(evidence));
        this.supplier = new SAMSupplier(minMappingQuality, bafRegions);
    }

    @NotNull
    public String contig() {
        return contig;
    }

    @NotNull
    public List<TumorBAF> evidence() {
        return evidence.stream().filter(x -> x.tumorIndelCount() == 0).collect(Collectors.toList());
    }

    @Override
    public TumorBAFEvidence call() throws Exception {

        try (SamReader reader = samReaderFactory.open(new File(bamFile))) {
            supplier.readOnce(reader, this::record);
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
