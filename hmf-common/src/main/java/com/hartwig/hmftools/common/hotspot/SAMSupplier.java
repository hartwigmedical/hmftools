package com.hartwig.hmftools.common.hotspot;

import java.util.Collection;
import java.util.Set;
import java.util.function.Consumer;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

class SAMSupplier {

    private static final int MIN_MAPPING_QUALITY = 1;

    private final ListMultimap<Chromosome, GenomeRegion> codingRegions;

    SAMSupplier(@NotNull final Collection<GenomeRegion> regions) {
        this.codingRegions = Multimaps.fromRegions(regions);
    }

    public void readOnce(@NotNull final SamReader samReader, @NotNull final Consumer<SAMRecord> consumer) {
        final Set<String> processed = Sets.newHashSet();
        for (GenomeRegion codingRegion : codingRegions.values()) {
            try (final SAMRecordIterator iterator = samReader.queryOverlapping(codingRegion.chromosome(),
                    (int) codingRegion.start(),
                    (int) codingRegion.end())) {
                while (iterator.hasNext()) {
                    final SAMRecord record = iterator.next();
                    if (samRecordMeetsQualityRequirements(record)) {
                        if (processed.add(record.toString())) {
                            consumer.accept(record);
                        }
                    }
                }
            }
        }
    }

    public boolean isInCodingRegions(@NotNull final GenomePosition hotspot) {
        return codingRegions.get(HumanChromosome.fromString(hotspot.chromosome())).stream().anyMatch(x -> x.contains(hotspot));
    }

    private boolean samRecordMeetsQualityRequirements(@NotNull final SAMRecord record) {
        return record.getMappingQuality() >= MIN_MAPPING_QUALITY && !record.getDuplicateReadFlag();
    }

}
