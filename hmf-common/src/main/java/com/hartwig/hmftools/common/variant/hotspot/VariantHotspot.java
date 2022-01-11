package com.hartwig.hmftools.common.variant.hotspot;

import static com.hartwig.hmftools.common.codon.Codons.isCodonMultiple;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantHotspot extends GenomePosition {

    @NotNull
    String ref();

    @NotNull
    String alt();

    default int end() {
        return position() + ref().length() - 1;
    }

    default boolean isSNV() {
        return ref().length() == 1 && alt().length() == 1;
    }

    default boolean isMNV() {
        return ref().length() == alt().length() && ref().length() != 1;
    }

    default boolean isIndel() {
        return ref().length() != alt().length();
    }

    default int indelLength()
    {
        if (ref().length() == alt().length())
            return 0;

        return alt().length() - ref().length();
    }
}
