package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextConsumerDispatcher implements Consumer<SAMRecord> {

    private final List<ReadContextConsumer> consumers;
    private final GenomePositionSelector<ReadContextConsumer> consumerSelector;

    public ReadContextConsumerDispatcher(@NotNull final List<ReadContextConsumer> consumerList) {
        consumerSelector = GenomePositionSelectorFactory.create(consumerList);
        this.consumers = consumerList;
    }

    @Override
    public void accept(final SAMRecord samRecord) {
        final GenomeRegion samRegion =
                GenomeRegions.create(samRecord.getContig(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());
        consumerSelector.select(samRegion, x -> x.accept(samRecord));
    }

    public List<ReadContextConsumer> consumers() {
        return consumers;
    }
}
