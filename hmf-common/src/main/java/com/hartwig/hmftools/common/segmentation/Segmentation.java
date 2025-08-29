package com.hartwig.hmftools.common.segmentation;

import java.util.Arrays;
import java.util.List;

public class Segmentation
{
    private final List<double[]> segments;

    public Segmentation(List<double[]> segments)
    {
        this.segments = segments;
    }

    public double[] completeSeries()
    {
        int length = 0;
        for(double[] segment : segments)
        {
            length += segment.length;
        }
        double[] result = new double[length];
        int index = 0;
        for(double[] segment : segments)
        {
            for(double value : segment)
            {
                result[index++] = value;
            }
        }
        return result;
    }

    public PCF pcf()
    {
        int cursor = 0;
        int[] lengths = new int[segments.size()];
        int[] starts = new int[segments.size()];
        double[] means = new double[segments.size()];
        for(int index = 0; index < segments.size(); index++)
        {
            double[] segment = segments.get(index);
            lengths[index] = segment.length;
            starts[index] = cursor;
            cursor += segment.length;
            means[index] = round(Stats.mean(segment));
        }
        return new PCF(lengths, starts, means);
    }

    public double cost(double gamma)
    {
        double sum = 0;
        for(double[] segment : segments)
        {
            sum += segmentCost(segment);
        }
        return sum + gamma * segments.size();
    }

    public double alternativeCost(double gamma)
    {
        double part1 = 0.0;
        double part2 = 0.0;
        for(double[] segment : segments)
        {
            double segmentSum = 0.0;
            for(double d : segment)
            {
                part1 += d * d;
                segmentSum += d;
            }
            part2 += (segmentSum * segmentSum) / segment.length;
        }
        return part1 - part2 + gamma * segments.size();
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder("Segmentation[");
        for(int i = 0; i < segments.size(); i++)
        {
            if(i > 0)
            {
                sb.append(",");
            }
            sb.append(Arrays.toString(segments.get(i)));
        }
        sb.append("]");
        return sb.toString();
    }

    private double segmentCost(double[] segment)
    {
        double mean = Stats.mean(segment);
        double sum = 0.0;
        for(double value : segment)
        {
            sum += Math.pow(value - mean, 2);
        }
        return sum;
    }

    private double round(double value)
    {
        return Math.round(value * 1000) / 1000.0;
    }

    @Override
    public boolean equals(Object o)
    {
        if(this == o)
        {
            return true;
        }
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        Segmentation that = (Segmentation) o;

        if(segments.size() != that.segments.size())
        {
            return false;
        }
        for(int i = 0; i < segments.size(); i++)
        {
            if(!Arrays.equals(segments.get(i), that.segments.get(i)))
            {
                return false;
            }
        }
        return true;
    }

    @Override
    public int hashCode()
    {
        int result = 17;
        for(double[] segment : segments)
        {
            result = result * 31 + Arrays.hashCode(segment);
        }
        return result;
    }
}