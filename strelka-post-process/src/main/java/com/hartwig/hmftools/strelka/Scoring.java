package com.hartwig.hmftools.strelka;

import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.strelka.scores.ImmutableReadScore;
import com.hartwig.hmftools.strelka.scores.ReadScore;
import com.hartwig.hmftools.strelka.scores.ReadType;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

final class Scoring {
    private static final Logger LOGGER = LogManager.getLogger(Scoring.class);

    private Scoring() {
    }

    private static ReadType getReadType(@NotNull final SAMRecord record, @NotNull final VariantContext variant) {
        // MIVO: assumes single alt allele
        final Allele alt = variant.getAlternateAllele(0);
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
            for (int index = 0; index < alt.length(); index++) {
                final int recordIndex = recordIdxOfVariantStart + index;
                if (record.getReadString().charAt(recordIndex - 1) == alt.getBaseString().charAt(index)) {
                    if (index != 0 && record.getReferencePositionAtReadPosition(recordIndex) != 0) {
                        return ReadType.REF;
                    }
                } else {
                    return ReadType.REF;
                }
            }
            // MIVO: check that this insertion is not part of a longer insertion
            if (record.getReadLength() >= recordIdxOfVariantStart + alt.length()
                    && record.getReferencePositionAtReadPosition(recordIdxOfVariantStart + alt.length()) == 0) {
                return ReadType.REF;
            }
            return ReadType.ALT;
        }
        if (variant.isSimpleDeletion()) {
            final Allele ref = variant.getReference();
            for (int index = 0; index < ref.length(); index++) {
                if (index != 0 && record.getReadPositionAtReferencePosition(variant.getStart() + index) != 0) {
                    return ReadType.REF;
                }
            }
            // MIVO: check that this deletion is not part of a longer deletion
            if (record.getAlignmentEnd() >= variant.getStart() + ref.length()
                    && record.getReadPositionAtReferencePosition(variant.getStart() + ref.length()) == 0) {
                return ReadType.REF;
            }
            return ReadType.ALT;
        }
        return ReadType.OTHER;
    }

    static ReadScore getReadScore(@NotNull final SAMRecord record, @NotNull final VariantContext variant) {
        final ReadType readType = getReadType(record, variant);
        // MIVO: assumes single alt allele
        final Allele alt = variant.getAlternateAllele(0);
        final int recordIdxOfVariantStart = record.getReadPositionAtReferencePosition(variant.getStart());
        switch (readType) {
            case REF: {
                if (variant.isSNP()) {
                    return ReadScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                } else if (variant.isSimpleInsertion()) {
                    return ReadScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                } else if (variant.isSimpleDeletion()) {
                    final int endIndex = Math.min(recordIdxOfVariantStart - 1 + variant.getReference().length(), record.getReadLength());
                    return ReadScore.of(readType, record.getBaseQualityString().substring(recordIdxOfVariantStart - 1, endIndex));
                }
                break;
            }
            case ALT: {
                if (variant.isSNP()) {
                    return ReadScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                } else if (variant.isSimpleInsertion()) {
                    return ReadScore.of(readType, record.getBaseQualityString()
                            .substring(recordIdxOfVariantStart - 1, recordIdxOfVariantStart - 1 + alt.length()));
                } else if (variant.isSimpleDeletion()) {
                    //MIVO: read score of next base after deletion if present, otherwise read score of base before deletion
                    if (record.getReadLength() > recordIdxOfVariantStart) {
                        return ReadScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart));
                    } else {
                        return ReadScore.of(readType, record.getBaseQualityString().charAt(recordIdxOfVariantStart - 1));
                    }
                }
                break;
            }
        }
        return ImmutableReadScore.of(readType, 0);
    }

    static Map<VariantContext, ReadScore> recordScores(@NotNull final SAMRecord record, @NotNull final PotentialMNV potentialMnv) {
        return potentialMnv.variants()
                .stream()
                .map(variant -> ImmutablePair.of(variant, getReadScore(record, variant)))
                .collect(Collectors.toMap(Pair::getLeft, Pair::getRight));
    }

    @NotNull
    @VisibleForTesting
    static Map<SAMRecord, Map<VariantContext, ReadScore>> correctMissingScores(
            @NotNull final Map<SAMRecord, Map<VariantContext, ReadScore>> scores) {
        final Map<VariantContext, Integer> averageScoresPerPosition = scores.entrySet()
                .stream()
                .flatMap(readEntry -> readEntry.getValue().entrySet().stream())
                .filter(entry -> refOrAlt(entry.getValue().type()))
                .collect(Collectors.groupingBy(Map.Entry::getKey, Collectors.mapping(readScore -> readScore.getValue().score(),
                        Collectors.collectingAndThen(Collectors.toList(),
                                list -> (int) (list.stream().mapToLong(value -> (long) value).sum() / list.size())))));
        return scores.entrySet()
                .stream()
                .map(readEntry -> correctEntry(readEntry, averageScoresPerPosition))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    @NotNull
    private static Map.Entry<SAMRecord, Map<VariantContext, ReadScore>> correctEntry(
            @NotNull final Map.Entry<SAMRecord, Map<VariantContext, ReadScore>> entry,
            @NotNull final Map<VariantContext, Integer> averageScoresPerPosition) {
        final Map.Entry<SAMRecord, Map<VariantContext, ReadScore>> correctedEntry = ImmutablePair.of(entry.getKey(), entry.getValue()
                .entrySet()
                .stream()
                .map(variantEntry -> ImmutablePair.of(variantEntry.getKey(),
                        correctRead(variantEntry.getValue(), averageScoresPerPosition.get(variantEntry.getKey()))))
                .filter(variantEntry -> refOrAlt(variantEntry.getRight().type()))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));
        final int readsWithOther = correctedEntry.getValue().entrySet().size() - entry.getValue().entrySet().size();
        if (readsWithOther > 0) {
            LOGGER.warn("Filtered out {} reads where read type was OTHER", readsWithOther);
        }
        return correctedEntry;
    }

    @NotNull
    private static ReadScore correctRead(@NotNull final ReadScore score, final int averageScoreAtPosition) {
        if (score.type() == ReadType.MISSING) {
            return ImmutableReadScore.of(ReadType.REF, averageScoreAtPosition);
        }
        return score;
    }

    private static boolean refOrAlt(@NotNull final ReadType readType) {
        return readType == ReadType.ALT || readType == ReadType.REF;
    }
}
