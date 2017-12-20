package com.hartwig.hmftools.common.variant;

import java.util.List;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(overshadowImplementation = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SomaticVariantImpl implements SomaticVariant {

    @Override
    @NotNull
    @Value.Default
    public VariantType type() {
        return VariantType.UNDEFINED;
    }

    @Override
    @NotNull
    @Value.Default
    public String chromosome() {
        return Strings.EMPTY;
    }

    @Override
    @Value.Default
    public long position() {
        return 0L;
    }

    @Override
    @NotNull
    @Value.Default
    public String ref() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    @Value.Default
    public String alt() {
        return Strings.EMPTY;
    }

    @Override
    @NotNull
    @Value.Default
    public String filter() {
        return Strings.EMPTY;
    }

    @Nullable
    @Value.Default
    public String dbsnpID() {
        return null;
    }

    @Nullable
    @Value.Default
    public String cosmicID() {
        return null;
    }

    @NotNull
    public abstract List<VariantAnnotation> annotations();

    @NotNull
    public abstract List<String> callers();

    @Override
    @Value.Default
    public int totalReadCount() {
        return 0;
    }

    @Override
    @Value.Default
    public int alleleReadCount() {
        return 0;
    }

    @Override
    @Value.Default
    public double mappability() {
        return 0;
    }

    @Override
    public String toString() {
        return "SomaticVariant{" + "chromosome='" + chromosome() + '\'' + ", position=" + position() + '}';
    }

    public static class Builder extends ImmutableSomaticVariantImpl.Builder implements VariantBuilder {

    }
}
