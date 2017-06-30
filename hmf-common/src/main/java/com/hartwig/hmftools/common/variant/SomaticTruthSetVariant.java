
package com.hartwig.hmftools.common.variant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SomaticTruthSetVariant implements Variant {

    @NotNull
    private final String originalVCFLine;
    
    @NotNull
    private final VariantType type;
    @NotNull
    private final String chromosome;
    private final long position;
    @NotNull
    private final String ref;
    @NotNull
    private final String alt;
    @NotNull
    private final String filter;

    public SomaticTruthSetVariant(@NotNull final String originalVCFLine, @NotNull final VariantType type, 
            @NotNull final String chromosome, final long position, @NotNull final String ref, 
            @NotNull final String alt, @NotNull final String filter) {
        this.originalVCFLine = originalVCFLine;
        this.type = type;
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.filter = filter;
    }
    
    @NotNull
    String originalVCFLine() {
        return originalVCFLine;
    }
    
    @Override
    public String ref() {
        return ref;
    }

    @Override
    public String alt() {
        return alt;
    }

    @Override
    public VariantType type() {
        return type;
    }

    @Override
    public String filter() {
        return filter;
    }

    @Override
    public String chromosome() {
        return chromosome;
    }

    @Override
    public long position() {
        return position;
    }
    
    @Override
    public String toString() {
        return "SomaticTruthSetVariant{" + "chromosome='" + chromosome + '\'' + ", position=" + position + '}';
    }
    
    public static class Builder implements VariantBuilder<SomaticTruthSetVariant> {
        @NotNull
        private String originalVCFLine = Strings.EMPTY;

        @NotNull
        private VariantType type = VariantType.UNDEFINED;
        @NotNull
        private String chromosome = Strings.EMPTY;
        private long position = 0L;
        @NotNull
        private String ref = Strings.EMPTY;
        @NotNull
        private String alt = Strings.EMPTY;
        @NotNull
        private String filter = Strings.EMPTY;
        
        @NotNull
        static SomaticTruthSetVariant.Builder fromVCF(@NotNull final String vcfLine) {
            final SomaticTruthSetVariant.Builder builder = new SomaticTruthSetVariant.Builder();
            builder.originalVCFLine(vcfLine);
            return builder;
        }

        @NotNull
        public static SomaticTruthSetVariant.Builder fromVariant(@NotNull final SomaticVariant variant) {
            final SomaticTruthSetVariant.Builder builder = fromVCF(variant.originalVCFLine());
            builder.type(variant.type());
            builder.chromosome(variant.chromosome());
            builder.position(variant.position());
            builder.ref(variant.ref());
            builder.alt(variant.alt());
            builder.filter(variant.filter());
            return builder;
        }

        public Builder() {
        }

        private void originalVCFLine(@NotNull final String originalVCFLine) {
            this.originalVCFLine = originalVCFLine;
        }

        @NotNull
        public SomaticTruthSetVariant.Builder type(@NotNull final VariantType type) {
            this.type = type;
            return this;
        }

        @NotNull
        public SomaticTruthSetVariant.Builder chromosome(@NotNull final String chromosome) {
            this.chromosome = chromosome;
            return this;
        }

        @NotNull
        public SomaticTruthSetVariant.Builder position(final long position) {
            this.position = position;
            return this;
        }

        @NotNull
        public SomaticTruthSetVariant.Builder ref(@NotNull final String ref) {
            this.ref = ref;
            return this;
        }

        @NotNull
        public SomaticTruthSetVariant.Builder alt(@NotNull final String alt) {
            this.alt = alt;
            return this;
        }

        @NotNull
        public SomaticTruthSetVariant.Builder filter(@NotNull final String filter) {
            this.filter = filter;
            return this;
        }
        
        @NotNull
        @Override
        public SomaticTruthSetVariant build() {
            return new SomaticTruthSetVariant(originalVCFLine, type, chromosome, position, ref, alt, filter);
        }
    }
}
