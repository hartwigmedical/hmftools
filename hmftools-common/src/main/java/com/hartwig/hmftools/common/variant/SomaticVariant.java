package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class SomaticVariant implements Variant {

    @NotNull
    private final VariantType type;
    @NotNull
    private final String filter;
    @NotNull
    private final String chromosome;
    private final long position;

    @NotNull
    private final List<String> callers;
    private final double alleleFrequency;
    private final boolean isDBSNP;
    private final boolean isCOSMIC;

    // KODU: Not sure if one variant can have multiple consequences.
    @NotNull
    private final List<VariantConsequence> consequences;

    private SomaticVariant(@NotNull final VariantType type, @NotNull final String filter,
            @NotNull final String chromosome, final long position, @NotNull final List<String> callers,
            final double alleleFrequency, final boolean isDBSNP, final boolean isCOSMIC,
            @NotNull final List<VariantConsequence> consequences) {
        this.type = type;
        this.filter = filter;
        this.chromosome = chromosome;
        this.position = position;
        this.callers = callers;
        this.alleleFrequency = alleleFrequency;
        this.isDBSNP = isDBSNP;
        this.isCOSMIC = isCOSMIC;
        this.consequences = consequences;
    }

    @NotNull
    public VariantType type() {
        return type;
    }

    @NotNull
    public String filter() {
        return filter;
    }

    @NotNull
    public String chromosome() {
        return chromosome;
    }

    public long position() {
        return position;
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

    public boolean isDBSNP() {
        return isDBSNP;
    }

    public boolean isCOSMIC() {
        return isCOSMIC;
    }

    public boolean isMissense() {
        return consequences.contains(VariantConsequence.MISSENSE_VARIANT);
    }

    @Override
    public String toString() {
        return "SomaticVariant{" + "chromosome='" + chromosome + '\'' + ", position=" + position + '}';
    }

    public static class Builder {
        @NotNull
        private final VariantType type;
        @NotNull
        private String filter = Strings.EMPTY;
        @NotNull
        private String chromosome = Strings.EMPTY;
        private long position = 0L;

        @NotNull
        private List<String> callers = Lists.newArrayList();
        private double alleleFrequency = Double.NaN;
        private boolean isDBSNP = false;
        private boolean isCOSMIC = false;
        @NotNull
        private List<VariantConsequence> consequences = Lists.newArrayList();

        public Builder(@NotNull final VariantType type) {
            this.type = type;
        }

        @NotNull
        public Builder filter(@NotNull final String filter) {
            this.filter = filter;
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
        public Builder isDBSNP(final boolean isDBSNP) {
            this.isDBSNP = isDBSNP;
            return this;
        }

        @NotNull
        public Builder isCOSMIC(final boolean isCOSMIC) {
            this.isCOSMIC = isCOSMIC;
            return this;
        }

        @NotNull
        Builder consequences(@NotNull final List<VariantConsequence> consequences) {
            this.consequences = consequences;
            return this;
        }

        @NotNull
        public SomaticVariant build() {
            return new SomaticVariant(type, filter, chromosome, position, callers, alleleFrequency, isDBSNP, isCOSMIC,
                    consequences);
        }
    }
}
