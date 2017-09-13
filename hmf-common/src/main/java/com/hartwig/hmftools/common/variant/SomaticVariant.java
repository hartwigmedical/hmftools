package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticVariant implements Variant, AllelicDepth {

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
    @Nullable
    private final String dbsnpID;
    @Nullable
    private final String cosmicID;
    @NotNull
    private final List<VariantAnnotation> annotations;
    @NotNull
    private final List<String> callers;
    private final int totalReadCount;
    private final int alleleReadCount;
    @NotNull
    private final String genoType;

    private SomaticVariant(@NotNull final String originalVCFLine, @NotNull final VariantType type,
            @NotNull final String chromosome, final long position, @NotNull final String ref,
            @NotNull final String alt, @NotNull final String filter, @Nullable final String dbsnpID,
            @Nullable final String cosmicID, @NotNull final List<VariantAnnotation> annotations,
            @NotNull final List<String> callers, final int totalReadCount, final int alleleReadCount,
            @NotNull String genoType) {
        this.originalVCFLine = originalVCFLine;
        this.type = type;
        this.chromosome = chromosome;
        this.position = position;
        this.ref = ref;
        this.alt = alt;
        this.filter = filter;
        this.dbsnpID = dbsnpID;
        this.cosmicID = cosmicID;
        this.annotations = annotations;
        this.callers = callers;
        this.totalReadCount = totalReadCount;
        this.alleleReadCount = alleleReadCount;
        this.genoType = genoType;
    }

    @NotNull
    String originalVCFLine() {
        return originalVCFLine;
    }

    @NotNull
    public VariantType type() {
        return type;
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    public long position() {
        return position;
    }

    @NotNull
    public String ref() {
        return ref;
    }

    @NotNull
    public String alt() {
        return alt;
    }

    @NotNull
    public String filter() {
        return filter;
    }

    @Nullable
    public String dbsnpID() {
        return dbsnpID;
    }

    public boolean isDBSNP() {
        return dbsnpID != null;
    }

    public boolean isCOSMIC() {
        return cosmicID != null;
    }

    @Nullable
    public String cosmicID() {
        return cosmicID;
    }

    @NotNull
    public List<VariantAnnotation> annotations() {
        return annotations;
    }

    @NotNull
    public String genoType() {
        return genoType;
    }
    
    public boolean hasConsequence(@NotNull final VariantConsequence consequence) {
        for (final VariantAnnotation annotation : annotations) {
            if (annotation.consequences().contains(consequence)) {
                return true;
            }
        }
        return false;
    }

    @NotNull
    public List<String> callers() {
        return callers;
    }

    public long callerCount() {
        return callers.size();
    }

    @Override
    public int totalReadCount() {
        return totalReadCount;
    }

    @Override
    public int alleleReadCount() {
        return alleleReadCount;
    }

    @Override
    public String toString() {
        return "SomaticVariant{" + "chromosome='" + chromosome + '\'' + ", position=" + position + '}';
    }

    public static class Builder implements VariantBuilder<SomaticVariant> {
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
        @Nullable
        private String dbsnpID = null;
        @Nullable
        private String cosmicID = null;
        @NotNull
        private List<VariantAnnotation> annotations = Lists.newArrayList();
        @NotNull
        private List<String> callers = Lists.newArrayList();
        private int totalReadCount = 0;
        private int alleleReadCount = 0;
        @NotNull
        private String genoType = Strings.EMPTY;

        @NotNull
        static Builder fromVCF(@NotNull final String vcfLine) {
            final Builder builder = new Builder();
            builder.originalVCFLine(vcfLine);
            return builder;
        }

        @NotNull
        public static Builder fromVariant(@NotNull final SomaticVariant variant) {
            final Builder builder = fromVCF(variant.originalVCFLine());
            builder.type(variant.type());
            builder.chromosome(variant.chromosome());
            builder.position(variant.position());
            builder.ref(variant.ref());
            builder.alt(variant.alt());
            builder.filter(variant.filter());
            builder.dnsnpID(variant.dbsnpID());
            builder.cosmicID(variant.cosmicID());
            builder.annotations(variant.annotations());
            builder.callers(variant.callers());
            builder.alleleReadCount(variant.alleleReadCount());
            builder.totalReadCount(variant.totalReadCount());
            builder.genoType(variant.genoType());
            return builder;
        }

        public Builder() {
        }

        private void originalVCFLine(@NotNull final String originalVCFLine) {
            this.originalVCFLine = originalVCFLine;
        }

        @NotNull
        public Builder type(@NotNull final VariantType type) {
            this.type = type;
            return this;
        }

        @NotNull
        public Builder chromosome(@NotNull final String chromosome) {
            this.chromosome = chromosome;
            return this;
        }

        @NotNull
        public Builder position(final long position) {
            this.position = position;
            return this;
        }

        @NotNull
        public Builder ref(@NotNull final String ref) {
            this.ref = ref;
            return this;
        }

        @NotNull
        public Builder alt(@NotNull final String alt) {
            this.alt = alt;
            return this;
        }

        @NotNull
        public Builder filter(@NotNull final String filter) {
            this.filter = filter;
            return this;
        }

        @NotNull
        public Builder dnsnpID(@Nullable final String dbsnpID) {
            this.dbsnpID = dbsnpID;
            return this;
        }

        @NotNull
        public Builder cosmicID(@Nullable final String cosmicID) {
            this.cosmicID = cosmicID;
            return this;
        }

        @NotNull
        public Builder annotations(@NotNull final List<VariantAnnotation> annotations) {
            this.annotations = annotations;
            return this;
        }

        @NotNull
        public Builder callers(@NotNull final List<String> callers) {
            this.callers = callers;
            return this;
        }

        @NotNull
        public Builder totalReadCount(final int totalReadCount) {
            this.totalReadCount = totalReadCount;
            return this;
        }

        @NotNull
        public Builder alleleReadCount(final int alleleReadCount) {
            this.alleleReadCount = alleleReadCount;
            return this;
        }
        
        public Builder genoType(@NotNull final String genoType) {
            this.genoType = genoType;
            return this;
        }

        @NotNull
        @Override
        public SomaticVariant build() {
            return new SomaticVariant(originalVCFLine, type, chromosome, position, ref, alt, filter, dbsnpID, cosmicID,
                    annotations, callers, totalReadCount, alleleReadCount, genoType);
        }
    }
}
