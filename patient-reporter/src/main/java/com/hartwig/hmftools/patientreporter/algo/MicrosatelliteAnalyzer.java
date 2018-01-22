package com.hartwig.hmftools.patientreporter.algo;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import com.hartwig.hmftools.common.purple.repeat.RepeatContext;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class },
             allParameters = true)
public abstract class MicrosatelliteAnalyzer {

    @NotNull
    public abstract IndexedFastaSequenceFile reference();

    @NotNull
    public static MicrosatelliteAnalyzer of(@NotNull final String fastaFileLocation) throws FileNotFoundException {
        return ImmutableMicrosatelliteAnalyzer.of(new IndexedFastaSequenceFile(new File(fastaFileLocation)));
    }

    public double analyzeVariants(@NotNull final List<SomaticVariant> variants) throws FileNotFoundException {
        double indelCount = 0;
        for (final SomaticVariant variant : variants) {
            if (isPassIndel(variant)) {
                long positionBeforeEvent = variant.position();
                long start = Math.max(positionBeforeEvent - 100, 1);
                long maxEnd = reference().getSequenceDictionary().getSequence(variant.chromosome()).getSequenceLength() - 1;
                long end = Math.min(positionBeforeEvent + 100, maxEnd);
                int relativePosition = (int) (positionBeforeEvent - start);
                final String sequence = reference().getSubsequenceAt(variant.chromosome(), start, end).getBaseString();
                if (EnrichedSomaticVariantFactory.getRepeatContext(variant, relativePosition, sequence)
                        .filter(this::repeatContextIsRelevant)
                        .isPresent()) {
                    indelCount++;
                }
            }
        }
        return indelCount / 3095;
    }

    private boolean isPassIndel(@NotNull final SomaticVariant variant) {
        return variant.filter().equals("PASS") && variant.ref().length() != variant.alt().length() && variant.ref().length() < 50
                && variant.alt().length() < 50;
    }

    private boolean repeatContextIsRelevant(@NotNull final RepeatContext repeatContext) {
        final int repeatCount = repeatContext.count();
        final int repeatSequenceLength = repeatContext.sequence().length();
        return repeatCount > 0 && ((repeatSequenceLength >= 2 && repeatSequenceLength <= 4) || (repeatSequenceLength == 1
                && repeatCount >= 5));
    }
}
