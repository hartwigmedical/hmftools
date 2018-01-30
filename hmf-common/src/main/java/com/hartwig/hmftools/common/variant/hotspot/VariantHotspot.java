package com.hartwig.hmftools.common.variant.hotspot;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantHotspot extends GenomePosition, Predicate<SomaticVariant> {

    @NotNull
    String ref();

    @NotNull
    String alt();

    @Override
    default boolean test(SomaticVariant variant) {
        return ref().equals(variant.ref()) && alt().equals(variant.alt()) && chromosome().equals(variant.chromosome())
                && position() == variant.position();
    }
}
