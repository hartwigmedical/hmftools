package com.hartwig.hmftools.common.segmentation.copynumber;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;

class Segmentation
{
    private final List<double[]> segments;

    Segmentation(List<double[]> segments)
    {
        this.segments = segments;
    }

    double[] completeSeries()
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

    PiecewiseConstantFit pcf()
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
            means[index] = round(Doubles.mean(segment));
        }
        return new PiecewiseConstantFit(lengths, starts, means);
    }

    double cost(double gamma)
    {
        double sum = 0;
        for(double[] segment : segments)
        {
            sum += segmentCost(segment);
        }
        return sum + gamma * segments.size();
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

    private double segmentCost(double[] segment)
    {
        if(segment.length == 0)
        {
            return 0.0;
        }

        // Calculate mean and sum of squares in a single pass
        double sum = 0.0;
        double sumSquared = 0.0;

        for(double value : segment)
        {
            sum += value;
            sumSquared += value * value;
        }

        double mean = sum / segment.length;

        // sum((x - mean)²) = sum(x²) - 2*mean*sum(x) + n*mean²
        // Simplified to: sum(x²) - n*mean²
        return sumSquared - segment.length * mean * mean;
    }

    private double round(double value)
    {
        return Math.round(value * 1000) / 1000.0;
    }
}