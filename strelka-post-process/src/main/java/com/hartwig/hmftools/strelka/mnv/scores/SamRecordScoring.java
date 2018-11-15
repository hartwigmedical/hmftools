package com.hartwig.hmftools.strelka.mnv.scores;

import java.util.Collection;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.sam.SamRecords;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public final class SamRecordScoring {
    private SamRecordScoring() {
    }

    @NotNull
    private static ReadType getReadType(@NotNull final SAMRecord record, @NotNull final VariantContext variant) {
        final Allele alt = variant.getAlternateAllele(0);
        final Allele ref = variant.getReference();
        final int recordIdxOfVariantStart = record.getReadPositionAtReferencePosition(variant.getStart());
        if (recordIdxOfVariantStart == 0) {
            //MIVO: variant position was deleted
            return ReadType.MISSING;

        }
        if (variant.isSNP()) {
            if (record.getReadString().charAt(recordIdxOfVariantStart - 1) == alt.getBaseString().charAt(0)) {
                return ReadType.ALT;
            } else {
                return ReadType.REF;
            }
        }
        if (variant.isSimpleInsertion()) {
            return SamRecords.containsInsert(variant.getStart(), alt.getBaseString(), record) ? ReadType.ALT : ReadType.REF;
        }
        if (variant.isSimpleDeletion()) {
            return SamRecords.containsDelete(variant.getStart(), ref.getBaseString(), record) ? ReadType.ALT : ReadType.REF;
        }
        return ReadType.OTHER;
    }

    // MIVO: assumes single alt allele
    @NotNull
    @VisibleForTesting
    static VariantScore getVariantScore(@NotNull final SAMRecord record, @NotNull final VariantContext variant) {
        assert variant.getAlternateAlleles().size() == 1;
        final ReadType readType = getReadType(record, variant);
        final Allele alt = variant.getAlternateAllele(0);
        final int recordIdxOfVariantStart = record.getReadPositionAtReferencePosition(variant.getStart());
        switch (readType) {
            case REF: {
                if (variant.isSNP() || variant.isSimpleInsertion()) {
                    return VariantScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                } else if (variant.isSimpleDeletion()) {
                    final int endIndex = Math.min(recordIdxOfVariantStart - 1 + variant.getReference().length(), record.getReadLength());
                    return VariantScore.of(readType, record.getBaseQualityString().substring(recordIdxOfVariantStart - 1, endIndex));
                }
                break;
            }
            case ALT: {
                if (variant.isSNP()) {
                    return VariantScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                } else if (variant.isSimpleInsertion()) {
                    return VariantScore.of(readType,
                            record.getBaseQualityString()
                                    .substring(recordIdxOfVariantStart - 1, recordIdxOfVariantStart - 1 + alt.length()));
                } else if (variant.isSimpleDeletion()) {
                    //MIVO: read score of next base after deletion if present, otherwise read score of base before deletion
                    if (record.getReadLength() > recordIdxOfVariantStart) {
                        return VariantScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart));
                    } else {
                        return VariantScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                    }
                }
                break;
            }
        }
        return ImmutableVariantScore.of(readType, 0);
    }

    @NotNull
    public static Map<VariantContext, VariantScore> scoresPerVariant(@NotNull final SAMRecord record,
            @NotNull final Collection<VariantContext> variants) {
        return variants.stream()
                .map(variant -> ImmutablePair.of(variant, getVariantScore(record, variant)))
                .collect(Collectors.toMap(Pair::getLeft, Pair::getRight));
    }
}