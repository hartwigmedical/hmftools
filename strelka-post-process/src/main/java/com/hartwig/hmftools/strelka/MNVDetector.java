package com.hartwig.hmftools.strelka;

import static com.hartwig.hmftools.strelka.VariantContextUtils.splitMultiAlleleVariant;

import java.io.File;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

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
    //TODO: MIVO: export potential mnv regions to .bed file

    private static final Logger LOGGER = LogManager.getLogger(MNVDetector.class);

    @NotNull
    protected abstract String tumorBAM();

    @Value.Derived
    @NotNull
    protected SamReader tumorReader() {
        return SamReaderFactory.makeDefault().open(new File(tumorBAM()));
    }

    @Value.Derived
    @NotNull
    protected SAMSequenceDictionary tumorDictionary() {
        return tumorReader().getFileHeader().getSequenceDictionary();
    }

    @NotNull
    Pair<PotentialMNVRegion, List<VariantContext>> checkMNV(@NotNull final PotentialMNVRegion potentialMnvRegion,
            @NotNull final VariantContext variant) {
        if (potentialMnvRegion.chromosome().equals(variant.getContig()) && variant.getStart() - potentialMnvRegion.end() <= 1
                && variant.getStart() - potentialMnvRegion.start() >= 0) {
            return ImmutablePair.of(PotentialMNVRegion.addVariant(potentialMnvRegion, variant), Lists.newArrayList());
        } else {
            return ImmutablePair.of(PotentialMNVRegion.fromVariant(variant), mergeVariants(potentialMnvRegion));
        }
    }

    @NotNull
    List<VariantContext> mergeVariants(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        if (potentialMnvRegion.potentialMnvs().size() == 0) {
            return potentialMnvRegion.variants();
        } else {
            LOGGER.info("Potential mnv of size {}: {}", potentialMnvRegion.variants().size(), potentialMnvRegion);
            final SAMRecordIterator samIterator = queryBam(potentialMnvRegion);
            final MNVRegionValidator regionValidator = validateMNVs(samIterator, potentialMnvRegion);
            samIterator.close();
            final List<VariantContext> mergedVariants = outputVariants(regionValidator);
            LOGGER.info("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -");
            return mergedVariants;
        }
    }

    @NotNull
    private SAMRecordIterator queryBam(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        final int referenceIndex = getReferenceIndex(potentialMnvRegion);
        final QueryInterval[] queryIntervals =
                new QueryInterval[] { new QueryInterval(referenceIndex, potentialMnvRegion.start(), potentialMnvRegion.end() - 1) };
        return tumorReader().queryOverlapping(queryIntervals);
    }

    @NotNull
    @VisibleForTesting
    static MNVRegionValidator validateMNVs(@NotNull final Iterator<SAMRecord> iterator,
            @NotNull final PotentialMNVRegion potentialMnvRegion) {
        MNVRegionValidator regionValidator = MNVRegionValidator.of(potentialMnvRegion);
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            regionValidator = regionValidator.addSamRecord(record);
        }
        return regionValidator;
    }

    @NotNull
    @VisibleForTesting
    static List<VariantContext> outputVariants(@NotNull final MNVRegionValidator regionValidator) {
        final List<VariantContext> result = Lists.newArrayList();
        LOGGER.info("Number of valid mnvs: {}", regionValidator.validMnvs().size());
        for (final PotentialMNV validMnv : regionValidator.validMnvs()) {
            final VariantContext mergedVariant = MNVMerger.mergeVariants(validMnv, regionValidator.mostFrequentReads());
            if (validMnv.variants().size() > 2) {
                LOGGER.info("Region with {} variants; mergedVariant: {}", validMnv.variants().size(), mergedVariant);
            }
            result.add(mergedVariant);
            LOGGER.info("Percentage of reads with all variants: {}", regionValidator.mnvScores().get(validMnv).frequency());
        }
        result.addAll(nonMnvVariants(regionValidator));
        result.sort(Comparator.comparing(VariantContext::getStart).thenComparing(variantContext -> variantContext.getReference().length()));
        return result;
    }

    @NotNull
    private static List<VariantContext> nonMnvVariants(@NotNull final MNVRegionValidator regionValidator) {
        return nonMnvVariants(regionValidator.region(), regionValidator.validMnvs());
    }

    @NotNull
    @VisibleForTesting
    static List<VariantContext> nonMnvVariants(@NotNull final PotentialMNVRegion region, @NotNull final Set<PotentialMNV> validMnvs) {
        final List<VariantContext> snvs = Lists.newArrayList();
        final Set<VariantContext> variantsInValidMnvs =
                validMnvs.stream().flatMap(mnv -> mnv.variants().stream()).collect(Collectors.toSet());
        region.variants().forEach(variant -> {
            if (variant.getAlternateAlleles().size() > 1) {
                final List<VariantContext> notContainedAlts = splitMultiAlleleVariant(variant).stream()
                        .filter(alt -> !variantsContain(variantsInValidMnvs, alt))
                        .collect(Collectors.toList());
                if (notContainedAlts.size() == variant.getAlternateAlleles().size()) {
                    snvs.add(variant);
                } else {
                    snvs.addAll(notContainedAlts);
                }
            } else if (!variantsContain(variantsInValidMnvs, variant)) {
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
}
