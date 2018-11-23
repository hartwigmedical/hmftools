package com.hartwig.hmftools.amber;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.amber.ModifiableTumorBAF;
import com.hartwig.hmftools.common.amber.NormalBAF;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.amber.TumorBAFFactory;
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

public class TumorEvidence implements Callable<TumorEvidence> {

    private final String contig;
    private final String bamFile;
    private final SamReaderFactory samReaderFactory;
    private final List<GenomeRegion> bafRegions;
    private final List<ModifiableTumorBAF> evidence;
    private final GenomePositionSelector<ModifiableTumorBAF> selector;

    public TumorEvidence(final String contig, final String bamFile, final SamReaderFactory samReaderFactory,
            final List<NormalBAF> bafRegions) {
        this.contig = contig;
        this.bamFile = bamFile;
        this.samReaderFactory = samReaderFactory;

        final GenomeRegionBuilder builder = new GenomeRegionBuilder(contig, 1);
        for (NormalBAF bafRegion : bafRegions) {
            builder.addPosition(bafRegion.position());
        }

        this.bafRegions = builder.build();
        this.evidence = bafRegions.stream().map(TumorBAFFactory::create).collect(Collectors.toList());
        this.selector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(evidence));
    }

    public String contig() {
        return contig;
    }

    public List<TumorBAF> evidence() {
        return new ArrayList<>(evidence);
    }

    @Override
    public TumorEvidence call() throws Exception {

        final SAMSupplier supplier = new SAMSupplier(bafRegions);
        try (SamReader reader = samReaderFactory.open(new File(bamFile))) {
            supplier.readOnce(reader, this::record);
        }

        return this;
    }

    private void record(@NotNull final SAMRecord record) {
        selector.select(asRegion(record), bafEvidence -> TumorBAFFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record) {
        return GenomeRegionFactory.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }
}
