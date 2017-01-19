package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticVariant implements Variant {

    @NotNull
    private final String originalVCFLine;

    @NotNull
    private final VariantType type;
    @NotNull
    private final String chromosome;
    private final long position;
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
    private final double alleleFrequency;
    private final int readDepth;

    public SomaticVariant(@NotNull final String originalVCFLine, @NotNull final VariantType type,
            @NotNull final String chromosome, final long position, @NotNull final String filter,
            @Nullable final String dbsnpID, @Nullable final String cosmicID,
            @NotNull final List<VariantAnnotation> annotations, @NotNull final List<String> callers,
            final double alleleFrequency, final int readDepth) {
        this.originalVCFLine = originalVCFLine;
        this.type = type;
        this.chromosome = chromosome;
        this.position = position;
        this.filter = filter;
        this.dbsnpID = dbsnpID;
        this.cosmicID = cosmicID;
        this.annotations = annotations;
        this.callers = callers;
        this.alleleFrequency = alleleFrequency;
        this.readDepth = readDepth;
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
    public String filter() {
        return filter;
    }

    public boolean isDBSNP() {
        return dbsnpID != null;
    }

    public boolean isCOSMIC() {
        return cosmicID != null;
    }

    @NotNull
    public List<VariantAnnotation> annotations() {
        return annotations;
    }

    public boolean hasConsequence(@NotNull VariantConsequence consequence) {
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

    public double alleleFrequency() {
        return alleleFrequency;
    }

    int readDepth() {
        return readDepth;
    }

    @Override
    public String toString() {
        return "SomaticVariant{" + "chromosome='" + chromosome + '\'' + ", position=" + position + '}';
    }

    public static class Builder {
        @NotNull
        private String originalVCFLine = Strings.EMPTY;

        @NotNull
        private final VariantType type;
        @NotNull
        private String chromosome = Strings.EMPTY;
        private long position = 0L;
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
        private double alleleFrequency = Double.NaN;
        private int readDepth = 0;

        @NotNull
        static Builder fromVCF(@NotNull final String vcfLine, @NotNull final VariantType type) {
            final Builder builder = new Builder(type);
            builder.originalVCFLine(vcfLine);
            return builder;
        }

        public Builder(@NotNull final VariantType type) {
            this.type = type;
        }

        private void originalVCFLine(@NotNull final String originalVCFLine) {
            this.originalVCFLine = originalVCFLine;
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
        Builder alleleFrequency(final double alleleFrequency) {
            this.alleleFrequency = alleleFrequency;
            return this;
        }

        @NotNull
        Builder readDepth(final int readDepth) {
            this.readDepth = readDepth;
            return this;
        }

        @NotNull
        public SomaticVariant build() {
            return new SomaticVariant(originalVCFLine, type, chromosome, position, filter, dbsnpID, cosmicID,
                    annotations, callers, alleleFrequency, readDepth);
        }
    }
}
