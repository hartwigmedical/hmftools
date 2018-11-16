package com.hartwig.hmftools.common.hotspot;

import java.util.Collection;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.base.Strings;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class InframeIndelHotspots {

    private final IndexedFastaSequenceFile sequenceFile;
    private final ListMultimap<Chromosome, GenomeRegion> codingRegions;

    public InframeIndelHotspots(@NotNull final Collection<GenomeRegion> regions, @NotNull final IndexedFastaSequenceFile sequenceFile) {
        this.sequenceFile = sequenceFile;
        this.codingRegions = Multimaps.fromRegions(regions);
    }

    @NotNull
    public Set<VariantHotspot> findInframeIndels(@NotNull final SamReader samReader) {
        final Set<VariantHotspot> indelsWithIncorrectRefs = Sets.newHashSet();

        for (GenomeRegion codingRegion : codingRegions.values()) {
            try (final SAMRecordIterator iterator = samReader.queryOverlapping(codingRegion.chromosome(),
                    (int) codingRegion.start(),
                    (int) codingRegion.end())) {
                while (iterator.hasNext()) {
                    final SAMRecord record = iterator.next();
                    indelsWithIncorrectRefs.addAll(findInframeIndelsWithIncorrectRefs(record));
                }
            }
        }

        return indelsWithIncorrectRefs.stream()
                .filter(this::isInCodingRegions)
                .map(this::correctRef)
                .filter(x -> x.isSimpleDelete() || x.isSimpleInsert())
                .collect(Collectors.toSet());
    }

    private boolean isInCodingRegions(@NotNull final VariantHotspot hotspot) {
        return codingRegions.get(HumanChromosome.fromString(hotspot.chromosome())).stream().anyMatch(x -> x.contains(hotspot));
    }

    @NotNull
    private VariantHotspot correctRef(@NotNull final VariantHotspot hotspot) {

        final ReferenceSequence refSequence = hotspot.isSimpleInsert()
                ? sequenceFile.getSubsequenceAt(hotspot.chromosome(), hotspot.position(), hotspot.position())
                : sequenceFile.getSubsequenceAt(hotspot.chromosome(), hotspot.position(), hotspot.position() + hotspot.ref().length() - 1);

        return ImmutableVariantHotspot.builder().from(hotspot).ref(refSequence.getBaseString()).build();
    }

    @NotNull
    static Set<VariantHotspot> findInframeIndelsWithIncorrectRefs(@NotNull final SAMRecord record) {
        Set<VariantHotspot> result = Sets.newHashSet();
        if (containsInframeIndel(record)) {
            for (int refPosition = record.getAlignmentStart(); refPosition <= record.getAlignmentEnd(); refPosition++) {
                int readPosition = record.getReadPositionAtReferencePosition(refPosition);

                if (readPosition != 0) {
                    int basesInserted = SAMRecords.basesInsertedAfterPosition(record, refPosition);
                    if (basesInserted > 0 && basesInserted % 3 == 0) {
                        result.add(createInsert(record, readPosition, refPosition, basesInserted + 1));
                        continue;
                    }

                    int basesDeleted = SAMRecords.basesDeletedAfterPosition(record, refPosition);
                    if (basesDeleted > 0 && basesDeleted % 3 == 0) {
                        result.add(createDelete(record, readPosition, refPosition, basesDeleted + 1));
                    }
                }
            }
        }

        return result;
    }

    @NotNull
    private static VariantHotspot createInsert(@NotNull final SAMRecord record, int readPosition, int refPosition, int length) {
        final String alt = record.getReadString().substring(readPosition - 1, readPosition - 1 + length);
        final String ref = alt.substring(0, 1);
        return ImmutableVariantHotspot.builder().chromosome(record.getContig()).position(refPosition).ref(ref).alt(alt).build();
    }

    @NotNull
    private static VariantHotspot createDelete(@NotNull final SAMRecord record, int readPosition, int refPosition, int length) {
        final String alt = record.getReadString().substring(readPosition - 1, readPosition);
        final String ref = alt + Strings.repeat("N", length - 1);
        return ImmutableVariantHotspot.builder().chromosome(record.getContig()).position(refPosition).ref(ref).alt(alt).build();
    }

    private static boolean containsInframeIndel(@NotNull final SAMRecord record) {
        for (final CigarElement cigarElement : record.getCigar().getCigarElements()) {
            if (cigarElement.getOperator().isIndel() && cigarElement.getLength() % 3 == 0) {
                return true;
            }
        }

        return false;
    }

}
