package com.hartwig.hmftools.common.segmentation.copynumber;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

public class SegmentationTestBase
{
    public double[] d(Number... values)
    {
        double[] result = new double[values.length];
        for(int i = 0; i < values.length; i++)
        {
            result[i] = values[i].doubleValue();
        }
        return result;
    }

    Segmentation segmentation(double[]... arrays)
    {
        return new Segmentation(java.util.Arrays.asList(arrays));
    }

    public File loadFileFromResources(String filename)
    {
        URL resource = ClassLoader.getSystemResource("segmentation/" + filename);
        if(resource == null)
        {
            throw new IllegalArgumentException("File " + filename + " not found in resources.");
        }
        try
        {
            return new File(resource.toURI());
        }
        catch(Exception e)
        {
            throw new RuntimeException("Error loading file: " + filename, e);
        }
    }

    public List<String[]> readResourcesFile(String filename)
    {
        List<String[]> result = new ArrayList<>();
        File file = loadFileFromResources(filename);

        try(BufferedReader reader = new BufferedReader(new FileReader(file)))
        {
            String line;
            while((line = reader.readLine()) != null)
            {
                String separator = filename.endsWith("csv") ? "," : "\t";
                result.add(line.split(separator));
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Error reading file: " + filename, e);
        }
        return result;
    }

    public static double round4(double value)
    {
        return Math.round(value * 10000) / 10000.0;
    }

    public static double round10(double value)
    {
        return Math.round(value * 10_000_000_000L) / 10_000_000_000.0;
    }
}
