package com.hartwig.hmftools.breakpointinspector.datamodel;

import java.util.Set;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedVariantContext {

    private static final int SHORT_VARIANT_LENGTH = 1000;

    @NotNull
    public abstract VariantContext variant();

    @NotNull
    public abstract HMFStructuralVariantType type();

    @NotNull
    public abstract Location location1();

    @NotNull
    public abstract Location location2();

    @NotNull
    public abstract Range uncertainty1();

    @NotNull
    public abstract Range uncertainty2();

    public abstract boolean imprecise();

    public abstract boolean isTranslocation();

    @NotNull
    public abstract Set<String> filters();

    @NotNull
    public abstract String insertSequence();

    @NotNull
    public abstract String homologySequence();

    public abstract int orientationBP1();

    public abstract int orientationBP2();

    @Value.Derived
    public boolean isShortVariant() {
        boolean shortDelete = type() == HMFStructuralVariantType.DEL;
        shortDelete &= location1().ReferenceIndex == location2().ReferenceIndex;
        shortDelete &= (location2().Position - location1().Position) < SHORT_VARIANT_LENGTH;

        boolean shortDuplicate = type() == HMFStructuralVariantType.DUP;
        shortDuplicate &= location1().ReferenceIndex == location2().ReferenceIndex;
        shortDuplicate &= (location2().Position - location1().Position) < SHORT_VARIANT_LENGTH;

        return shortDelete || shortDuplicate;
    }

    @Value.Derived
    public boolean isInsert() {
        return type() == HMFStructuralVariantType.INS;
    }
}

