package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.Scoring.recordScores;
import static com.hartwig.hmftools.strelka.VariantContextUtils.splitMultiAlleleVariant;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
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

    Pair<PotentialMNVRegion, List<VariantContext>> checkMNV(@NotNull final PotentialMNVRegion potentialMnvRegion,
            @NotNull final VariantContext variant) {
        if (potentialMnvRegion.chromosome().equals(variant.getContig()) && variant.getStart() - potentialMnvRegion.end() <= 1
                && variant.getStart() - potentialMnvRegion.start() >= 0) {
            return ImmutablePair.of(PotentialMNVRegion.addVariant(potentialMnvRegion, variant), Lists.newArrayList());
        } else {
            return ImmutablePair.of(PotentialMNVRegion.fromVariant(variant), mergeVariants(potentialMnvRegion));
        }
    }

    List<VariantContext> mergeVariants(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        if (potentialMnvRegion.potentialMnvs().size() == 0) {
            return potentialMnvRegion.variants();
        } else {
            LOGGER.info("Potential mnv of size {}: {}", potentialMnvRegion.variants().size(), potentialMnvRegion);
            final List<VariantContext> mergedVariants = queryBam(potentialMnvRegion);
            LOGGER.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
            return mergedVariants;
        }
    }

    private List<VariantContext> queryBam(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        final List<VariantContext> result = Lists.newArrayList();
        final int referenceIndex = getReferenceIndex(potentialMnvRegion);
        final QueryInterval[] queryIntervals =
                new QueryInterval[] { new QueryInterval(referenceIndex, potentialMnvRegion.start(), potentialMnvRegion.end() - 1) };
        final SAMRecordIterator iterator = tumorReader().queryOverlapping(queryIntervals);
        final Map<Integer, GapReads> readsPerPosition =
                potentialMnvRegion.gapPositions().stream().collect(Collectors.toMap(Function.identity(), (position) -> GapReads.empty()));
        Map<PotentialMNV, MNVScore> scores = potentialMnvRegion.potentialMnvs()
                .stream()
                .collect(Collectors.toMap(Function.identity(), potentialMNV -> MNVScore.of(potentialMNV.variants())));
        //        MNVScore score = MNVScore.of(potentialMnvRegion.variants());
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            if (goodRead(record) && containsAllMNVPositions(record, potentialMnvRegion)) {
                final Map<VariantContext, ReadScore> samRecordScores = recordScores(record, potentialMnvRegion.variants());
                scores = scores.entrySet()
                        .stream()
                        .map(entry -> Pair.of(entry.getKey(), MNVScore.addReadScore(entry.getValue(), samRecordScores)))
                        .collect(Collectors.toMap(Pair::getKey, Pair::getValue));
                //score = MNVScore.addReadScore(score, samRecordScores);
                potentialMnvRegion.gapPositions()
                        .forEach(position -> readsPerPosition.put(position,
                                GapReads.addRead(readsPerPosition.get(position), getReadAtPosition(record, position))));
            }
        }

        if (potentialMnvRegion.variants().size() == 2) {
            final Map<Integer, Character> mostFrequentReadPerPosition = readsPerPosition.entrySet()
                    .stream()
                    .map(entry -> ImmutablePair.of(entry.getKey(), entry.getValue().mostFrequentRead()))
                    .collect(Collectors.toMap(Pair::getKey, Pair::getValue));
            for (final Map.Entry<PotentialMNV, MNVScore> scoreEntry : scores.entrySet()) {
                final MNVScore score = scoreEntry.getValue();
                LOGGER.info(readsPerPosition);
                if (score.isMNV()) {
                    //TODO: pass sample name as param?
                    final String sampleName = potentialMnvRegion.variants().get(0).getSampleNamesOrderedByName().get(0);
                    final VariantContext mergedVariant =
                            MNVMerger.mergeVariants(sampleName, potentialMnvRegion.variants(), mostFrequentReadPerPosition);
                    result.add(mergedVariant);
                } else {
                    result.addAll(potentialMnvRegion.variants());
                }
                LOGGER.info("Percentage of reads with both variants: {}", score.frequency());
            }
        }
        iterator.close();
        return result;
    }

    static List<VariantContext> nonMnvVariants(@NotNull final PotentialMNVRegion mnvRegion,
            @NotNull final Map<PotentialMNV, MNVScore> mnvScores) {
        final List<VariantContext> snvs = Lists.newArrayList();
        final Set<VariantContext> variantsInMnvs = mnvScores.entrySet()
                .stream()
                .filter(mnvEntry -> mnvEntry.getValue().isMNV())
                .flatMap(mnvEntry -> mnvEntry.getKey().variants().stream())
                .collect(Collectors.toSet());
        mnvRegion.variants().forEach(variant -> {
            if (variant.getAlternateAlleles().size() > 1) {
                final List<VariantContext> notContainedAlts = splitMultiAlleleVariant(variant).stream()
                        .filter(alt -> !variantsContain(variantsInMnvs, alt))
                        .collect(Collectors.toList());
                if (notContainedAlts.size() == variant.getAlternateAlleles().size()) {
                    snvs.add(variant);
                } else {
                    snvs.addAll(notContainedAlts);
                }
            } else if (!variantsContain(variantsInMnvs, variant)) {
                snvs.add(variant);
            }
        });
        return snvs;
    }

    private static boolean variantsContain(@NotNull final Collection<VariantContext> variants, @NotNull final VariantContext other) {
        for (final VariantContext variant : variants) {
            if (variant.getContig().equals(other.getContig()) && variant.getStart() == other.getStart()
                    && variant.getEnd() == other.getEnd() && variant.getReference().equals(other.getReference())
                    && variant.getAlternateAlleles().stream().allMatch(other.getAlternateAlleles()::contains)) {
                return true;
            }
        }
        return false;
    }

    private boolean goodRead(@NotNull final SAMRecord record) {
        return !record.getDuplicateReadFlag();
    }

    private static boolean containsAllMNVPositions(@NotNull final SAMRecord record, @NotNull final PotentialMNVRegion potentialMnvRegion) {
        final VariantContext lastVariant = potentialMnvRegion.variants().get(potentialMnvRegion.variants().size() - 1);
        return record.getAlignmentStart() <= potentialMnvRegion.start()
                && record.getAlignmentEnd() >= lastVariant.getStart() + lastVariant.getReference().length() - 1;
    }

    private int getReferenceIndex(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        int referenceIndex = tumorDictionary().getSequenceIndex(potentialMnvRegion.chromosome());
        if (referenceIndex >= 0) {
            return referenceIndex;
        }
        if (!potentialMnvRegion.chromosome().startsWith("chr")) {
            referenceIndex = tumorDictionary().getSequenceIndex("chr" + potentialMnvRegion.chromosome());
        } else {
            referenceIndex = tumorDictionary().getSequenceIndex(potentialMnvRegion.chromosome().substring(3));
        }
        if (referenceIndex < 0) {
            throw new RuntimeException(potentialMnvRegion.chromosome() + " is not in the BAM: " + tumorBAM());
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
