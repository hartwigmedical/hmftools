package com.hartwig.hmftools.common.hotspot;

import static com.hartwig.hmftools.common.sam.SAMRecords.basesDeletedAfterPosition;
import static com.hartwig.hmftools.common.sam.SAMRecords.basesInsertedAfterPosition;
import static com.hartwig.hmftools.common.sam.SAMRecords.containsDelete;
import static com.hartwig.hmftools.common.sam.SAMRecords.containsInsert;
import static com.hartwig.hmftools.common.sam.SAMRecords.getAvgBaseQuality;
import static com.hartwig.hmftools.common.sam.SAMRecords.getBaseQuality;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.sam.SAMRecords;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class VariantHotspotEvidenceFactory {

    static final int MIN_BASE_QUALITY = 13;
    private static final int SNV_MNV_BUFFER = 2;

    private final IndexedFastaSequenceFile sequenceFile;
    private final SamReader samReader;
    private final SAMSupplier samSupplier;

    public VariantHotspotEvidenceFactory(@NotNull final Collection<GenomeRegion> regions,
            @NotNull final IndexedFastaSequenceFile sequenceFile, @NotNull final SamReader samReader) {
        this.sequenceFile = sequenceFile;
        this.samReader = samReader;
        samSupplier = new SAMSupplier(regions);
    }

    @NotNull
    public List<VariantHotspotEvidence> evidence(@NotNull final Collection<VariantHotspot> hotspots) {

        final Map<VariantHotspot, String> refSequenceMap = Maps.newHashMap();
        final Map<VariantHotspot, ModifiableVariantHotspotEvidence> evidenceMap = Maps.newHashMap();
        final ListMultimap<Chromosome, VariantHotspot> hotspotMap = Multimaps.fromPositions(hotspots);

        final Consumer<SAMRecord> samRecordConsumer = record -> {
            final Chromosome samChromosome = HumanChromosome.fromString(record.getContig());

            if (hotspotMap.containsKey(samChromosome)) {

                for (VariantHotspot hotspot : hotspotMap.get(samChromosome)) {
                    int length = hotspot.ref().length();
                    int hotspotStartPosition = (int) hotspot.position();
                    int hotspotEndPosition = (int) hotspot.position() + length - 1;

                    if (samRecordOverlapsVariant(hotspotStartPosition, hotspotEndPosition, record) && startPositionValid(hotspot, record)) {
                        final ModifiableVariantHotspotEvidence evidence =
                                evidenceMap.computeIfAbsent(hotspot, VariantHotspotEvidenceFactory::create);
                        if (hotspot.isSimpleDelete()) {
                            findEvidenceOfDelete(evidence, hotspot, record);
                        } else if (hotspot.isSimpleInsert()) {
                            findEvidenceOfInsert(evidence, hotspot, record);
                        } else if (hotspot.isMNV() || hotspot.isSNV()) {
                            int sequenceLength = sequenceFile.getSequenceDictionary().getSequence(hotspot.chromosome()).getSequenceLength();
                            int bufferStartPosition = Math.max(0, hotspotStartPosition - SNV_MNV_BUFFER);
                            int bufferEndPosition = Math.min(sequenceLength, hotspotEndPosition + SNV_MNV_BUFFER);
                            if (samRecordOverlapsVariant(bufferStartPosition, bufferEndPosition, record)) {
                                final String refSequence = refSequenceMap.computeIfAbsent(hotspot,
                                        x -> refSequence(bufferStartPosition, bufferEndPosition, record.getContig()));
                                findEvidenceOfMNV(evidence, bufferStartPosition, refSequence, hotspot, record);
                            }
                        }
                    }
                }
            }
        };

        samSupplier.readOnce(samReader, samRecordConsumer);
        return new ArrayList<>(evidenceMap.values());
    }

    @NotNull
    private String refSequence(int start, int end, String contig) {
        return sequenceFile.getSubsequenceAt(contig, start, end).getBaseString();
    }

    @NotNull
    static ModifiableVariantHotspotEvidence findEvidenceOfInsert(@NotNull final ModifiableVariantHotspotEvidence builder,
            @NotNull final VariantHotspot hotspot, @NotNull final SAMRecord record) {
        assert (hotspot.isSimpleInsert());

        int hotspotStartPosition = (int) hotspot.position();
        int recordStartPosition = record.getReadPositionAtReferencePosition(hotspotStartPosition);
        if (recordStartPosition == 0) {
            return builder;
        }

        int recordStartQuality = getBaseQuality(record, recordStartPosition);
        if (containsInsert(record, hotspotStartPosition, hotspot.alt())) {
            int quality = getAvgBaseQuality(record, recordStartPosition, hotspot.alt().length());
            if (quality < MIN_BASE_QUALITY) {
                return builder;
            }
            return builder.setReadDepth(builder.readDepth() + 1)
                    .setIndelSupport(builder.indelSupport() + 1)
                    .setAltQuality(builder.altQuality() + quality)
                    .setAltSupport(builder.altSupport() + 1);
        }

        if (recordStartQuality < MIN_BASE_QUALITY) {
            return builder;
        }

        int insertedBases = basesInsertedAfterPosition(record, hotspotStartPosition);
        int deletedBases = basesDeletedAfterPosition(record, hotspotStartPosition);
        if (insertedBases == 0 && deletedBases == 0 && record.getReadString().charAt(recordStartPosition - 1) == hotspot.ref().charAt(0)) {
            builder.setRefSupport(builder.refSupport() + 1);
        } else if (insertedBases > 0 || deletedBases > 0) {
            builder.setIndelSupport(builder.indelSupport() + 1);
        }

        return builder.setReadDepth(builder.readDepth() + 1);
    }

    static ModifiableVariantHotspotEvidence findEvidenceOfDelete(@NotNull final ModifiableVariantHotspotEvidence builder,
            @NotNull final VariantHotspot hotspot, @NotNull final SAMRecord record) {
        assert (hotspot.isSimpleDelete());

        int hotspotStartPosition = (int) hotspot.position();
        int recordStartPosition = record.getReadPositionAtReferencePosition(hotspotStartPosition);
        if (recordStartPosition == 0) {
            return builder;
        }

        int recordStartQuality = getBaseQuality(record, recordStartPosition);

        if (containsDelete(record, hotspotStartPosition, hotspot.ref())) {
            int quality =
                    record.getReadLength() > recordStartPosition ? getAvgBaseQuality(record, recordStartPosition, 2) : recordStartQuality;
            if (quality < MIN_BASE_QUALITY) {
                return builder;
            }
            return builder.setReadDepth(builder.readDepth() + 1)
                    .setIndelSupport(builder.indelSupport() + 1)
                    .setAltQuality(builder.altQuality() + quality)
                    .setAltSupport(builder.altSupport() + 1);
        }

        if (recordStartQuality < MIN_BASE_QUALITY) {
            return builder;
        }

        int insertedBases = basesInsertedAfterPosition(record, hotspotStartPosition);
        int deletedBases = basesDeletedAfterPosition(record, hotspotStartPosition);
        if (insertedBases == 0 && deletedBases == 0 && record.getReadString().charAt(recordStartPosition - 1) == hotspot.ref().charAt(0)) {
            builder.setRefSupport(builder.refSupport() + 1);
        } else if (insertedBases > 0 || deletedBases > 0) {
            builder.setIndelSupport(builder.indelSupport() + 1);
        }

        return builder.setReadDepth(builder.readDepth() + 1);
    }

    @NotNull
    static ModifiableVariantHotspotEvidence findEvidenceOfMNV(@NotNull final ModifiableVariantHotspotEvidence builder, int start,
            @NotNull final String refSequence, @NotNull final VariantHotspot hotspot, @NotNull final SAMRecord record) {

        int hotspotStartPosition = (int) hotspot.position();
        int hotspotLength = Math.max(hotspot.ref().length(), hotspot.alt().length());

        int recordStartPosition = record.getReadPositionAtReferencePosition(hotspotStartPosition);
        int recordStartQuality = SAMRecords.getBaseQuality(record, recordStartPosition);

        if (isVariantPartOfLargerMNV(start, refSequence, hotspot, record)) {
            return recordStartQuality < MIN_BASE_QUALITY ? builder : builder.setReadDepth(builder.readDepth() + 1);
        }

        for (int i = 0; i < hotspotLength; i++) {
            int readPosition = record.getReadPositionAtReferencePosition(hotspotStartPosition + i);
            boolean isDeleted = readPosition == 0;
            boolean isInserted = record.getReferencePositionAtReadPosition(recordStartPosition + i) == 0;

            if (isInserted || isDeleted) {
                if (recordStartQuality < MIN_BASE_QUALITY) {
                    return builder;
                }
                return builder.setReadDepth(builder.readDepth() + 1);
            }
        }

        final String samBases = record.getReadString().substring(recordStartPosition - 1, recordStartPosition - 1 + hotspotLength);
        if (samBases.equals(hotspot.alt())) {
            int altQuality = getAvgBaseQuality(record, recordStartPosition, hotspotLength);
            if (altQuality < MIN_BASE_QUALITY) {
                return builder;
            }

            builder.setAltQuality(builder.altQuality() + altQuality).setAltSupport(builder.altSupport() + 1);
        } else if (samBases.equals(hotspot.ref())) {
            builder.setRefSupport(builder.refSupport() + 1);
        }

        return builder.setReadDepth(builder.readDepth() + 1);
    }

    private boolean startPositionValid(@NotNull final VariantHotspot hotspot, @NotNull final SAMRecord record) {
        return record.getReadPositionAtReferencePosition((int) hotspot.position()) != 0;
    }

    private boolean samRecordOverlapsVariant(int start, int end, @NotNull final SAMRecord record) {
        return record.getAlignmentStart() <= start && record.getAlignmentEnd() >= end;
    }

    @VisibleForTesting
    static boolean isVariantPartOfLargerMNV(int start, @NotNull final String refSequence, @NotNull final VariantHotspot hotspot,
            @NotNull final SAMRecord record) {

        int mvnLength = Math.max(hotspot.ref().length(), hotspot.alt().length());
        int hotspotStartPosition = (int) hotspot.position();
        int hotspotEndPosition = (int) hotspot.position() + mvnLength - 1;

        int requiredBuffer = hotspotStartPosition - start;
        if (requiredBuffer == 0) {
            return false;
        }

        return isStartPartOfLargerMNV(start, record, refSequence.substring(0, requiredBuffer)) || isEndPartOfLargerMNV(
                hotspotEndPosition + 1, record, refSequence.substring(requiredBuffer + mvnLength));

    }

    private static boolean isStartPartOfLargerMNV(int samOffset, @NotNull final SAMRecord samRecord, @NotNull final String refSequence) {
        int variantStart = samOffset + refSequence.length();

        int startReadPosition = samRecord.getReadPositionAtReferencePosition(variantStart);
        assert (startReadPosition != 0);

        for (int i = 1; i <= refSequence.length(); i++) {
            int bufferReadPosition = samRecord.getReadPositionAtReferencePosition(variantStart - i);

            if (bufferReadPosition == 0) {
                return false; // DEL
            }

            if (samRecord.getReferencePositionAtReadPosition(startReadPosition - i) == 0) {
                return false; // INS
            }

            int refSequenceIndex = refSequence.length() - i;
            if (samRecord.getReadString().charAt(bufferReadPosition - 1) != refSequence.charAt(refSequenceIndex)) {
                return true;
            }

        }

        return false;
    }

    private static boolean isEndPartOfLargerMNV(int samOffset, @NotNull final SAMRecord samRecord, @NotNull final String refSequence) {
        int variantEnd = samOffset - 1;

        int endReadPosition = samRecord.getReadPositionAtReferencePosition(variantEnd);
        assert (endReadPosition != 0);

        for (int i = 1; i <= refSequence.length(); i++) {
            int bufferReadPosition = samRecord.getReadPositionAtReferencePosition(variantEnd + i);

            if (bufferReadPosition == 0) {
                return false; // DEL
            }

            if (samRecord.getReferencePositionAtReadPosition(endReadPosition + i) == 0) {
                return false; // INS
            }

            if (samRecord.getReadString().charAt(bufferReadPosition - 1) != refSequence.charAt(i - 1)) {
                return true;
            }

        }

        return false;
    }

    @NotNull
    static ModifiableVariantHotspotEvidence create(@NotNull final VariantHotspot hotspot) {
        return ModifiableVariantHotspotEvidence.create()
                .from(hotspot)
                .setRef(hotspot.ref())
                .setAlt(hotspot.alt())
                .setAltQuality(0)
                .setAltSupport(0)
                .setRefSupport(0)
                .setIndelSupport(0)
                .setReadDepth(0);
    }

}
