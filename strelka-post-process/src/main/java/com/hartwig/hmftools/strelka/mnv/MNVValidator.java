package com.hartwig.hmftools.strelka.mnv;

import java.io.File;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

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
public abstract class MNVValidator {
    @NotNull
    abstract String tumorBAM();

    @NotNull
    @Value.Derived
    protected SamReader tumorReader() {
        return SamReaderFactory.makeDefault().open(new File(tumorBAM()));
    }

    @NotNull
    @Value.Derived
    protected SAMSequenceDictionary tumorDictionary() {
        return tumorReader().getFileHeader().getSequenceDictionary();
    }

    @NotNull
    public List<VariantContext> mergeVariants(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        if (potentialMnvRegion.potentialMnvs().size() == 0) {
            return potentialMnvRegion.variants();
        } else {
            final SAMRecordIterator samIterator = queryBam(potentialMnvRegion);
            final MNVRegionValidator regionValidator = validateMNVs(samIterator, potentialMnvRegion);
            samIterator.close();
            return outputVariants(regionValidator);
        }
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
        for (final PotentialMNV validMnv : regionValidator.validMnvs()) {
            final VariantContext mergedVariant = MNVMerger.mergeVariants(validMnv, regionValidator.mostFrequentReads());
            result.add(mergedVariant);
        }
        result.addAll(regionValidator.nonMnvVariants());
        result.sort(Comparator.comparing(VariantContext::getStart).thenComparing(variantContext -> variantContext.getReference().length()));
        return result;
    }

    @NotNull
    private SAMRecordIterator queryBam(@NotNull final PotentialMNVRegion potentialMnvRegion) {
        final int referenceIndex = getReferenceIndex(potentialMnvRegion);
        final QueryInterval[] queryIntervals =
                new QueryInterval[] { new QueryInterval(referenceIndex, potentialMnvRegion.start(), potentialMnvRegion.end() - 1) };
        return tumorReader().queryOverlapping(queryIntervals);
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
