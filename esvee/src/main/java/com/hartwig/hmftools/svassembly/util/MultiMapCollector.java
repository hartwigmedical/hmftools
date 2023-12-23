package com.hartwig.hmftools.svassembly.util;

import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;

public class MultiMapCollector<INPUT, KEY, VALUE> implements Collector<INPUT, Map<KEY, List<VALUE>>, Map<KEY, List<VALUE>>>
{
    private final Function<INPUT, KEY> mKeyExtractor;
    private final Function<INPUT, VALUE> mValueExtractor;

    public MultiMapCollector(final Function<INPUT, KEY> keyExtractor, final Function<INPUT, VALUE> valueExtractor)
    {
        mKeyExtractor = keyExtractor;
        mValueExtractor = valueExtractor;
    }

    public static <INPUT, KEY> MultiMapCollector<INPUT, KEY, INPUT> keyed(final Function<INPUT, KEY> keyExtractor)
    {
        return new MultiMapCollector<>(keyExtractor, Function.identity());
    }

    @Override
    public Supplier<Map<KEY, List<VALUE>>> supplier()
    {
        return HashMap::new;
    }

    @Override
    public BiConsumer<Map<KEY, List<VALUE>>, INPUT> accumulator()
    {
        return (map, input) -> map.computeIfAbsent(mKeyExtractor.apply(input), key -> new ArrayList<>()).add(mValueExtractor.apply(input));
    }

    @Override
    public BinaryOperator<Map<KEY, List<VALUE>>> combiner()
    {
        return (left, right) ->
        {
            right.forEach((key, value) -> left.merge(key, value, (l, r) ->
            {
                l.addAll(r);
                return l;
            }));
            return left;
        };
    }

    @Override
    public Function<Map<KEY, List<VALUE>>, Map<KEY, List<VALUE>>> finisher()
    {
        return Function.identity();
    }

    @Override
    public Set<Characteristics> characteristics()
    {
        return EnumSet.of(Characteristics.IDENTITY_FINISH, Characteristics.UNORDERED);
    }
}
