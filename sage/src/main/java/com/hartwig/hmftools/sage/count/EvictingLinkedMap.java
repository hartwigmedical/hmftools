package com.hartwig.hmftools.sage.count;

import java.util.AbstractMap;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Function;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EvictingLinkedMap<K, V> {

    private static final int MAX_CAPACITY = 1024;
    private static final int ACTUAL_CAPACITY = 700;

    private final Map<K, V> map;
    private final BiConsumer<K, V> evictionHandler;
    private final ArrayDeque<Map.Entry<K, V>> arrayDeque;

    public EvictingLinkedMap(final BiConsumer<K, V> evictionHandler) {
        this.evictionHandler = evictionHandler;
        this.map = new HashMap<>();
        arrayDeque = new ArrayDeque<>(MAX_CAPACITY);
    }

    @Nullable
    public V get(@NotNull final K key) {
        return map.get(key);
    }

    public V compute(K key, BiFunction<K, V, V> remapping) {

        BiFunction<K, V, V> internalRemapping = (k, oldValue) -> {

            V newValue = remapping.apply(k, oldValue);
            if (oldValue == null) {
                arrayDeque.addLast(new AbstractMap.SimpleEntry<>(key, newValue));

                if (arrayDeque.size() > ACTUAL_CAPACITY) {
                    evict(1);
                }
            }

            return newValue;
        };

        return map.compute(key, internalRemapping);
    }

    public V computeIfAbsent(K key, Function<K, V> supplier) {
        Function<K, V> internalSupplier = k -> {
            V newValue = supplier.apply(k);
            arrayDeque.addLast(new AbstractMap.SimpleEntry<>(key, newValue));
            if (arrayDeque.size() > ACTUAL_CAPACITY) {
                evict(1);
            }
            return newValue;
        };

        return map.computeIfAbsent(key, internalSupplier);
    }

    public void evictAll() {
        evict(arrayDeque.size());
    }

    private void evict(int count) {
        assert (count <= arrayDeque.size());

        for (int i = 1; i <= count; i++) {
            Map.Entry<K, V> entry = arrayDeque.removeFirst();
            map.remove(entry.getKey());
            evictionHandler.accept(entry.getKey(), entry.getValue());
        }
    }
}
