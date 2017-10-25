package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.Scoring.recordScores;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.strelka.scores.ReadScore;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class MNVDetector {
    //TODO: MIVO: deal with potential mnvs of 3+ variants where only a subset are actually mnvs
    //TODO: MIVO: multiple variants on same position => multiple potential mnvs starting at same position
    //TODO: MIVO: variants with multiple alts => multiple potential mnvs
    //TODO: MIVO: export potential mnv regions to .bed file

    private static final Logger LOGGER = LogManager.getLogger(MNVDetector.class);

    @NotNull
    protected abstract String tumorBAM();

    @Value.Derived
    protected SamReader tumorReader() {
        return SamReaderFactory.makeDefault().open(new File(tumorBAM()));
    }

    @Value.Derived
    protected SAMSequenceDictionary tumorDictionary() {
        return tumorReader().getFileHeader().getSequenceDictionary();
    }

    Pair<PotentialMNV, List<VariantContext>> checkMNV(@NotNull final PotentialMNV potentialMnv, @NotNull final VariantContext variant) {
        if (potentialMnv.chromosome().equals(variant.getContig()) && variant.getStart() - potentialMnv.end() <= 1
                && variant.getStart() - potentialMnv.end() >= 0) {
            return ImmutablePair.of(PotentialMNV.addVariant(potentialMnv, variant), Lists.newArrayList());
        } else {
            return ImmutablePair.of(PotentialMNV.fromVariant(variant), mergeVariants(potentialMnv));
        }
    }

    List<VariantContext> mergeVariants(@NotNull final PotentialMNV potentialMnv) {
        if (potentialMnv.variants().size() <= 1) {
            return potentialMnv.variants();
        } else {
            LOGGER.info("Potential mnv of size {}: {}", potentialMnv.variants().size(), potentialMnv);
            final List<VariantContext> mergedVariants = queryBam(potentialMnv);
            LOGGER.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
            return mergedVariants;
        }
    }

    private List<VariantContext> queryBam(@NotNull final PotentialMNV potentialMnv) {
        final List<VariantContext> result = Lists.newArrayList();
        final int referenceIndex = getReferenceIndex(potentialMnv);
        final QueryInterval[] queryIntervals =
                new QueryInterval[] { new QueryInterval(referenceIndex, potentialMnv.start(), potentialMnv.lastPosition()) };
        final SAMRecordIterator iterator = tumorReader().queryOverlapping(queryIntervals);
        final Map<Integer, GapReads> readsPerPosition =
                potentialMnv.gapPositions().stream().collect(Collectors.toMap(Function.identity(), (position) -> GapReads.empty()));
        MNVScore score = MNVScore.of(potentialMnv.variants());
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            if (goodRead(record) && containsAllMNVPositions(record, potentialMnv)) {
                final Map<VariantContext, ReadScore> samRecordScores = recordScores(record, potentialMnv);
                score = MNVScore.addReadScore(score, samRecordScores);
                potentialMnv.gapPositions()
                        .forEach(position -> readsPerPosition.put(position,
                                GapReads.addRead(readsPerPosition.get(position), getReadAtPosition(record, position))));
            }
        }
        if (potentialMnv.variants().size() == 2) {
            final Map<Integer, Character> mostFrequentReadPerPosition = readsPerPosition.entrySet()
                    .stream()
                    .map(entry -> ImmutablePair.of(entry.getKey(), entry.getValue().mostFrequentRead()))
                    .collect(Collectors.toMap(Pair::getKey, Pair::getValue));
            if (score.isMNV()) {
                //TODO: pass sample name as param?
                final String sampleName = potentialMnv.variants().get(0).getSampleNamesOrderedByName().get(0);
                final VariantContext mergedVariant =
                        MNVMerger.mergeVariants(sampleName, potentialMnv.variants(), mostFrequentReadPerPosition);
                result.add(mergedVariant);
            } else {
                result.addAll(potentialMnv.variants());
            }
            LOGGER.info("Percentage of reads with both variants: {}", score.frequency());
        }
        iterator.close();
        return result;
    }

    private boolean goodRead(@NotNull final SAMRecord record) {
        return !record.getDuplicateReadFlag();
    }

    private static boolean containsAllMNVPositions(@NotNull final SAMRecord record, @NotNull final PotentialMNV potentialMnv) {
        final VariantContext lastVariant = potentialMnv.variants().get(potentialMnv.variants().size() - 1);
        return record.getAlignmentStart() <= potentialMnv.start()
                && record.getAlignmentEnd() >= lastVariant.getStart() + lastVariant.getReference().length() - 1;
    }

    private int getReferenceIndex(@NotNull final PotentialMNV potentialMnv) {
        int referenceIndex = tumorDictionary().getSequenceIndex(potentialMnv.chromosome());
        if (referenceIndex >= 0) {
            return referenceIndex;
        }
        if (!potentialMnv.chromosome().startsWith("chr")) {
            referenceIndex = tumorDictionary().getSequenceIndex("chr" + potentialMnv.chromosome());
        } else {
            referenceIndex = tumorDictionary().getSequenceIndex(potentialMnv.chromosome().substring(3));
        }
        if (referenceIndex < 0) {
            throw new RuntimeException(potentialMnv.chromosome() + " is not in the BAM: " + tumorBAM());
        }
        return referenceIndex;
    }

    private static Character getReadAtPosition(@NotNull final SAMRecord record, final int position) {
        final int recordPosition = record.getReadPositionAtReferencePosition(position);
        if (recordPosition == 0) {
            return 'N';
        } else {
            return record.getReadString().charAt(recordPosition - 1);
        }
    }
}
