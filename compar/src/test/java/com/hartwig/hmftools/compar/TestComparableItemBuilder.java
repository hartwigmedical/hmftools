package com.hartwig.hmftools.compar;

import java.util.function.Consumer;
import java.util.function.Function;

public class TestComparableItemBuilder<B, T>
{
    private final SupplierWithException<B> builderSupplier;
    private final Function<B, T> buildFunction;
    private final Consumer<B> alternateDefaults;

    public TestComparableItemBuilder(final SupplierWithException<B> builderSupplier, final Function<B, T> buildFunction,
            final Consumer<B> alternateDefaults)
    {
        this.builderSupplier = builderSupplier;
        this.buildFunction = buildFunction;
        this.alternateDefaults = alternateDefaults;
    }

    public T create(final Consumer<B> initializer)
    {
        try
        {
            B builder = builderSupplier.get();
            initializer.accept(builder);
            return buildFunction.apply(builder);
        }
        catch(Exception e)
        {
            throw new RuntimeException("Failed to create test data object", e);
        }
    }

    public T create()
    {
        return create(b -> {});
    }

    public T createWithAlternateDefaults(final Consumer<B> initializer)
    {
        try
        {
            B builder = builderSupplier.get();
            alternateDefaults.accept(builder);
            initializer.accept(builder);
            return buildFunction.apply(builder);
        }
        catch(Exception e)
        {
            throw new RuntimeException("Failed to create test data object with alternate defaults", e);
        }
    }

    public T createWithAlternateDefaults()
    {
        return createWithAlternateDefaults(b -> {});
    }

    // Allows checked exception in builder supplier
    @FunctionalInterface
    public interface SupplierWithException<T>
    {
        T get() throws Exception;
    }
}
