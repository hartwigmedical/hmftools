package com.hartwig.hmftools.common.collect;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.collect.PartitionUtils.mergePartitions;

import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.function.BinaryOperator;
import java.util.function.Supplier;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.Nullable;

public class TwoDGridMap<T>
{
    private static class Point implements Comparable<Point>
    {
        public final int X;
        public final int Y;

        public Point(int x, int y)
        {
            X = x;
            Y = y;
        }

        @Override
        public int compareTo(final Point o)
        {
            if(X != o.X)
                return X - o.X;

            return Y - o.Y;
        }

        @Override
        public boolean equals(final Object o)
        {
            if(this == o)
                return true;

            if(!(o instanceof Point))
                return false;

            final Point point = (Point) o;
            return X == point.X && Y == point.Y;
        }

        @Override
        public int hashCode()
        {
            return X + 31 * Y;
        }
    }

    private final static BinaryOperator<List<Point>> LIST_ACCUMULATOR = (acc, x) ->
    {
        acc.addAll(x);
        return acc;
    };

    private final SortedMap<Point, T> mContainer;
    private final OneDGridMap<List<Point>> mKeysSortedByX;
    private final OneDGridMap<List<Point>> mKeysSortedByY;

    public TwoDGridMap()
    {
        mContainer = Maps.newTreeMap();
        mKeysSortedByX = new OneDGridMap<>();
        mKeysSortedByY = new OneDGridMap<>();
    }

    @Nullable
    public T get(int x, int y)
    {
        Point p = new Point(x, y);
        return mContainer.get(p);
    }

    @Nullable
    private T get(final Point key)
    {
        return mContainer.get(key);
    }

    public void put(int x, int y, T value)
    {
        Point p = new Point(x, y);
        put(p, value);
    }

    private void put(final Point key, T value)
    {
        boolean keyAlreadyExists = mContainer.containsKey(key);
        mContainer.put(key, value);

        if(keyAlreadyExists)
            return;

        List<Point> xKeys = mKeysSortedByX.get(key.X);
        if(xKeys == null)
        {
            xKeys = Lists.newArrayList();
            mKeysSortedByX.put(key.X, xKeys);
        }

        xKeys.add(key);

        List<Point> yKeys = mKeysSortedByY.get(key.Y);
        if(yKeys == null)
        {
            yKeys = Lists.newArrayList();
            mKeysSortedByY.put(key.Y, yKeys);
        }

        yKeys.add(key);
    }

    public void merge(int x, int y, final T value, final BinaryOperator<T> acc)
    {
        Point p = new Point(x, y);
        T currentValue = mContainer.get(p);
        if(currentValue == null)
        {
            put(p, value);
            return;
        }

        T newValue = acc.apply(currentValue, value);
        mContainer.put(p, newValue);
    }

    public List<TwoDGridMap<T>> partitionByXDist(int maxDistance)
    {
        List<List<Point>> keysPartitions = mKeysSortedByX.mergeValuesByDistance(maxDistance, LIST_ACCUMULATOR, Lists::newArrayList);
        return partitionFromKeyPartition(keysPartitions);
    }

    public List<TwoDGridMap<T>> partitionByYDist(int maxDistance)
    {
        List<List<Point>> keysPartitions = mKeysSortedByY.mergeValuesByDistance(maxDistance, LIST_ACCUMULATOR, Lists::newArrayList);
        return partitionFromKeyPartition(keysPartitions);
    }

    private List<TwoDGridMap<T>> partitionFromKeyPartition(final List<List<Point>> keyPartitions)
    {
        List<TwoDGridMap<T>> partitions = Lists.newArrayList();
        for(List<Point> keysPartition : keyPartitions)
        {
            TwoDGridMap<T> partition = new TwoDGridMap<>();
            for(Point key : keysPartition)
            {
                T value = mContainer.get(key);
                partition.put(key, value);
            }

            partitions.add(partition);
        }

        return partitions;
    }

    public List<TwoDGridMap<T>> partitionByXDistThenYDist(int maxDistance)
    {
        List<TwoDGridMap<T>> partitionsByX = partitionByXDist(maxDistance);
        List<TwoDGridMap<T>> finalPartitions = Lists.newArrayList();
        for(TwoDGridMap<T> xPartition : partitionsByX)
            finalPartitions.addAll(xPartition.partitionByYDist(maxDistance));

        return finalPartitions;
    }

    public List<List<T>> partitionValuesByDistance(int maxDistance)
    {
        List<List<T>> partitions = Lists.newArrayList();
        for(TwoDGridMap<T> rectPartition : partitionByXDistThenYDist(maxDistance))
        {
            List<Point> keys = Lists.newArrayList(rectPartition.mContainer.keySet());
            UnionFind<Integer> merger = new UnionFind<>();
            for(int i = 0; i < keys.size(); i++)
                merger.add(i);

            for(int i = 0; i < keys.size() - 1; i++)
            {
                Point key1 = keys.get(i);
                for(int j = i + 1; j < keys.size(); j++)
                {
                    Point key2 = keys.get(j);
                    int dist = abs(key1.X - key2.X) + abs(key1.Y - key2.Y);
                    if(dist <= maxDistance)
                        merger.merge(i, j);
                }
            }

            for(Set<Integer> indexPartition : merger.getPartitions())
            {
                List<T> partition = indexPartition.stream().map(i -> rectPartition.get(keys.get(i))).collect(Collectors.toList());
                partitions.add(partition);
            }
        }

        return partitions;
    }

    public List<T> mergeValuesByDistance(int maxDistance, final BinaryOperator<T> accumulator)
    {
        return mergeValuesByDistance(maxDistance, accumulator, null);
    }

    public List<T> mergeValuesByDistance(int maxDistance, final BinaryOperator<T> accumulator, @Nullable final Supplier<T> identity)
    {
        return mergePartitions(partitionValuesByDistance(maxDistance), accumulator, identity);
    }
}
