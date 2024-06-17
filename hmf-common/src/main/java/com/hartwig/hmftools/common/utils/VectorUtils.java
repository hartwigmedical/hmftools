package com.hartwig.hmftools.common.utils;

import static com.hartwig.hmftools.common.sigs.DataUtils.doublesEqual;

import java.util.List;

import com.google.common.collect.Lists;

public final class VectorUtils
{
    public static double sumVector(double[] vec)
    {
        double total = 0;
        for (double v : vec)
        {
            total += v;
        }

        return total;
    }

    public static void sumVectors(final double[] source, double[] dest)
    {
        if(source.length != dest.length)
            return;

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] += source[i];
        }
    }

    public static double[] vectorMultiply(final double[] vec1, final double[] vec2)
    {
        if(vec1.length != vec2.length)
            return null;

        double[] output = new double[vec1.length];
        for(int i = 0; i < vec1.length; ++i)
        {
            output[i] = vec1[i] * vec2[i];
        }

        return output;
    }

    public static void clear(double[] vec)
    {
        initVector(vec, 0);
    }

    public static void initVector(double[] vec, double value)
    {
        for(int i = 0; i < vec.length; ++i)
        {
            vec[i] = value;
        }
    }

    public static void vectorMultiply(double[] vec, double value)
    {
        for(int i = 0; i < vec.length; ++i)
        {
            vec[i] *= value;
        }
    }

    public static void copyVector(final double[] source, double[] dest)
    {
        if(source.length != dest.length)
            return;

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] = source[i];
        }
    }

    public static void addVector(final double[] source, double[] dest)
    {
        if(source.length != dest.length)
            return;

        for(int i = 0; i < source.length; ++i)
        {
            dest[i] += source[i];
        }
    }

    public static List<Integer> getSortedVectorIndices(final double[] data, boolean ascending)
    {
        // returns a list of indices into the original vector, being the sorted data list
        List<Integer> sortedList = Lists.newArrayList();

        for(int i = 0; i < data.length; ++i)
        {
            if(i == 0)
            {
                sortedList.add(i);
                continue;
            }

            int j = 0;
            for(; j < sortedList.size(); ++j)
            {
                int origIndex = sortedList.get(j);

                if(ascending && data[i] < data[origIndex])
                    break;
                else if(!ascending && data[i] > data[origIndex])
                    break;
            }

            sortedList.add(j, i);
        }

        return sortedList;
    }

    public static void optimisedAdd(final List<Double> items, final double value, boolean ascending)
    {
        // build a sorted list using a binary search

        // early exits
        if(items.isEmpty())
        {
            items.add(value);
            return;
        }

        if((ascending && value < items.get(0) || (!ascending && value > items.get(0))))
        {
            items.add(0, value);
            return;
        }

        int itemCount = items.size();
        if((ascending && value > items.get(itemCount - 1)) || (!ascending && value < items.get(itemCount - 1)))
        {
            items.add(value);
            return;
        }

        if(itemCount < 20)
        {
            int index = 0;
            while(index < items.size())
            {
                if(ascending && value < items.get(index))
                    break;
                if(!ascending && value > items.get(index))
                    break;

                ++index;
            }

            items.add(index, value);
            return;
        }

        int lowIndex = 0;
        int highIndex = items.size() - 1;
        int currentIndex = items.size() / 2;

        while(true)
        {
            double currentValue = items.get(currentIndex);

            if(currentValue == value)
            {
                items.add(currentIndex, value);
                return;
            }

            if((ascending && value < currentValue) || (!ascending && value > currentValue))
            {
                // current index is looking too high in the list
                if(currentIndex == lowIndex + 1)
                {
                    // no need to look any lower (again
                    items.add(currentIndex, value);
                    return;
                }

                highIndex = currentIndex;
            }
            else
            {
                if(currentIndex == highIndex - 1)
                {
                    items.add(currentIndex + 1, value);
                    return;
                }

                lowIndex = currentIndex;
            }

            int newIndex = lowIndex + (highIndex - lowIndex) / 2;

            if(newIndex == currentIndex)
            {
                items.add(currentIndex, value);
                return;
            }

            currentIndex = newIndex;
        }
    }
}
