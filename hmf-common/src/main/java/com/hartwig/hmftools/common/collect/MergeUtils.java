package com.hartwig.hmftools.common.collect;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiPredicate;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.commons.lang3.tuple.Pair;

public class MergeUtils
{
    public static <K, V> List<V> clusterMerger(final Map<K, V> elements, final BiPredicate<K, K> canMergeFn,
            final Comparator<V> valueComparator, final BinaryOperator<V> mergeFn)
    {
        if(elements.isEmpty())
            return Lists.newArrayList();

        if(elements.size() == 1)
            return Lists.newArrayList(elements.values());

        Comparator<Map.Entry<K, V>> heapComparator = (final Map.Entry<K, V> x, final Map.Entry<K, V> y) ->
                valueComparator.compare(x.getValue(), y.getValue());

        List<K> allKeys = Lists.newArrayList(elements.keySet());
        Map<K, Set<K>> keyAdjacency = Maps.newHashMap();
        for(K key : allKeys)
            keyAdjacency.put(key, Sets.newHashSet());

        for(int i = 0; i < allKeys.size() - 1; i++)
        {
            K key1 = allKeys.get(i);
            for(int j = i + 1; j < allKeys.size(); j++)
            {
                K key2 = allKeys.get(j);
                if(canMergeFn.test(key1, key2))
                {
                    keyAdjacency.get(key1).add(key2);
                    keyAdjacency.get(key2).add(key1);
                }
            }
        }

        List<V> mergedValues = Lists.newArrayList();
        Set<K> unprocessedKeys = Sets.newHashSet(allKeys);

        Heap<Map.Entry<K, V>> entryHeap = new Heap<>(heapComparator);
        entryHeap.addAll(elements.entrySet());
        while(!entryHeap.isEmpty())
        {
            Map.Entry<K, V> entry = entryHeap.pop();
            if(!unprocessedKeys.contains(entry.getKey()))
                continue;

            unprocessedKeys.remove(entry.getKey());
            V mergedValue = entry.getValue();
            Set<K> neighbours = keyAdjacency.get(entry.getKey()).stream()
                    .filter(unprocessedKeys::contains)
                    .collect(Collectors.toCollection(Sets::newHashSet));

            if(neighbours.isEmpty())
            {
                mergedValues.add(mergedValue);
                continue;
            }

            Heap<Map.Entry<K, V>> neighbourHeap = new Heap<>(heapComparator);
            for(K neighbour : neighbours)
            {
                Map.Entry<K, V> neighbourEntry = Pair.of(neighbour, elements.get(neighbour));
                neighbourHeap.add(neighbourEntry);
            }

            while(!neighbourHeap.isEmpty())
            {
                Map.Entry<K, V> neighbourEntry = neighbourHeap.pop();
                if(!neighbours.contains(neighbourEntry.getKey()))
                    continue;

                unprocessedKeys.remove(neighbourEntry.getKey());
                neighbours.retainAll(keyAdjacency.get(neighbourEntry.getKey()));
                mergedValue = mergeFn.apply(mergedValue, neighbourEntry.getValue());
            }

            mergedValues.add(mergedValue);
        }

        return mergedValues;
    }
}
