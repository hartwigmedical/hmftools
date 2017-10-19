package com.hartwig.hmftools.common.purple.segment;

import java.util.List;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class StructuralVariantCluster implements GenomeRegion {

    @NotNull
    public abstract List<StructuralVariantPosition> variants();

    public StructuralVariantSupport type() {
        if (variants().isEmpty()) {
            return StructuralVariantSupport.NONE;
        }

        if (variants().size() == 1) {
            return StructuralVariantSupport.fromVariant(variants().get(0).type());
        }

        return StructuralVariantSupport.MULTIPLE;
    }

    @Override
    public String toString() {

        StringBuilder builder = new StringBuilder();
        builder.append("Chromosome:")
                .append(chromosome())
                .append(", Start:")
                .append(start())
                .append(", End:")
                .append(end())
                .append(", Type:")
                .append(type())
                .append(", ")
                .append(variants().size())
                .append(" Variant(s):");
        for (StructuralVariantPosition sv : variants()) {
            builder.append(" [SV pos:")
                    .append(sv.position())
                    .append(", type:")
                    .append(sv.type())
                    .append(", orientation:")
                    .append(sv.orientation())
                    .append(", id:")
                    .append(sv.id())
                    .append("]");
        }

        return builder.toString();
    }

}
