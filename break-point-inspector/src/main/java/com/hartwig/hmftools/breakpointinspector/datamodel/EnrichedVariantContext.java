package com.hartwig.hmftools.breakpointinspector.datamodel;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnrichedVariantContext {

    public abstract VariantContext variant();

    public abstract HMFVariantType type();

    public abstract Location location1();

    public abstract Location location2();
}
