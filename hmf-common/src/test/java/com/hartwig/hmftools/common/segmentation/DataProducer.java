package com.hartwig.hmftools.common.segmentation;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

/**
 * Utility class for generating R code for testing.
 */
public class DataProducer
{

    @Test
    public void list1()
    {
        printRCode(1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    @Test
    public void list2()
    {
        printRCode(1.0, 1.0, 1.0, 20.0, 20.0, 20.0);
    }

    private void printRCode(double... doubles)
    {
        List<String> chromosomesList = new ArrayList<>();
        List<Integer> positionsList = new ArrayList<>();
        List<Double> doublesList = new ArrayList<>();
        int position = 1000000;
        for(double d : doubles)
        {
            chromosomesList.add("1");
            positionsList.add(position);
            position += 10000;
            doublesList.add(d);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("df = data.frame(\n");
        sb.append("  chrom = c(").append(String.join(", ", chromosomesList)).append("),\n");
        sb.append("  pos = c(").append(joinIntegers(positionsList)).append("),\n");
        sb.append("  blah = c(").append(joinDoubles(doublesList)).append(")\n");
        sb.append(")\n");
        sb.append("r1 <- pcf(df, kmin=1, gamma=0.0, verbose=TRUE, fast= FALSE)\n");
        sb.append("r1");
        System.out.println(sb.toString());

        StringBuilder sb2 = new StringBuilder();
        sb2.append("bl <- aaaexactPcf(c(").append(joinDoubles(doublesList)).append("), gamma=2.0, yest=TRUE)\n");
        sb2.append("bl");
        System.out.println(sb2.toString());
    }

    /**
     * Joins a list of integers into a comma-separated string.
     *
     * @param list The list to join
     * @return A comma-separated string
     */
    private String joinIntegers(List<Integer> list)
    {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < list.size(); i++)
        {
            if(i > 0)
            {
                sb.append(", ");
            }
            sb.append(list.get(i));
        }
        return sb.toString();
    }

    /**
     * Joins a list of doubles into a comma-separated string.
     *
     * @param list The list to join
     * @return A comma-separated string
     */
    private String joinDoubles(List<Double> list)
    {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < list.size(); i++)
        {
            if(i > 0)
            {
                sb.append(", ");
            }
            sb.append(list.get(i));
        }
        return sb.toString();
    }
}
