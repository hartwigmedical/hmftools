package com.hartwig.hmftools.common.collect;

import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.function.BiFunction;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class Cluster
{
    public static <T> int clusterCount(final List<T> data, BiFunction<T, T, Double> distanceFn, int maxDistance)
    {
        SortedMap<Integer, SortedSet<Integer>> neighbours = Maps.newTreeMap();

        for(int i = 0; i < data.size(); i++)
        {
            T data1 = data.get(i);
            neighbours.computeIfAbsent(i, key -> Sets.newTreeSet());
            for(int j = i + 1; j < data.size(); j++)
            {
                T data2 = data.get(j);
                neighbours.computeIfAbsent(j, key -> Sets.newTreeSet());
                if(distanceFn.apply(data1, data2) <= maxDistance)
                {
                    neighbours.get(i).add(j);
                    neighbours.get(j).add(i);
                }
            }
        }

        int count = 0;
        List<Integer> rootIndices = Lists.newArrayList(neighbours.keySet());
        for(int rootIdx : rootIndices)
        {
            if(neighbours.get(rootIdx).isEmpty())
            {
                neighbours.remove(rootIdx);
                count++;
            }
        }

        while(!neighbours.isEmpty())
        {
            int minNeighbourCount = Integer.MAX_VALUE;
            int minNeighbourIdx = 0;
            for(int rootIdx : neighbours.keySet())
            {
                int neighbourCount = neighbours.get(rootIdx).size();
                if(neighbourCount < minNeighbourCount)
                {
                    minNeighbourCount = neighbourCount;
                    minNeighbourIdx = rootIdx;
                }
            }

            count++;

            SortedSet<Integer> minNeighbours = neighbours.get(minNeighbourIdx);
            neighbours.remove(minNeighbourIdx);
            rootIndices = Lists.newArrayList(neighbours.keySet());
            for(int rootIdx : rootIndices)
            {
                if(minNeighbours.contains(rootIdx))
                {
                    neighbours.remove(rootIdx);
                    continue;
                }

                neighbours.get(rootIdx).remove(minNeighbourIdx);
                neighbours.get(rootIdx).removeAll(minNeighbours);
                if(neighbours.get(rootIdx).isEmpty())
                {
                    neighbours.remove(rootIdx);
                    count++;
                }
            }
        }

        return count;
    }
}
