package com.hartwig.hmftools.common.segmentation;

import org.junit.Test;
import static org.junit.Assert.assertEquals;

import java.util.List;

public class GammaTest extends SegmentationTestBase
{
    
    @Test
    public void compareWithR() {
        assertEquals(0.2361, gammaForFile("ratios1.tsv"), 0.0001);
        assertEquals(0.3325, gammaForFile("ratios4.tsv"), 0.0001);
        assertEquals(0.0, gammaForFile("ratios5.csv"), 0.0001);
        assertEquals(0.3099, gammaForFile("ratios_chr13_1.tsv"), 0.0001);
    }

    @Test
    public void filterWidthTest() {
        assertEquals(-1, Gamma.filterWidth(0));
        assertEquals(1, Gamma.filterWidth(1));
        assertEquals(1, Gamma.filterWidth(2));
        assertEquals(3, Gamma.filterWidth(3));
        assertEquals(3, Gamma.filterWidth(4));
        assertEquals(5, Gamma.filterWidth(5));
        assertEquals(5, Gamma.filterWidth(6));
        assertEquals(49, Gamma.filterWidth(49));
        assertEquals(49, Gamma.filterWidth(50));
        assertEquals(51, Gamma.filterWidth(51));
        assertEquals(51, Gamma.filterWidth(52));
        assertEquals(51, Gamma.filterWidth(52));
        assertEquals(51, Gamma.filterWidth(53));
        assertEquals(51, Gamma.filterWidth(54));
        assertEquals(51, Gamma.filterWidth(10000));
    }

    /**
     * Calculates the gamma value for a file.
     *
     * @param filename The name of the file
     * @param gamma The gamma parameter (default 50.0)
     * @param normalise Whether to normalise the segment penalty (default true)
     * @param column The column to use (default "ratio")
     * @return The gamma value
     */
    private double gammaForFile(String filename, double gamma, boolean normalise, String column) {
        List<String[]> data = readResourcesFile(filename);
        
        // Find the column index
        int columnIndex = -1;
        if (data.size() > 0) {
            String[] header = data.get(0);
            for (int i = 0; i < header.length; i++) {
                if (header[i].equals(column)) {
                    columnIndex = i;
                    break;
                }
            }
        }
        
        if (columnIndex == -1) {
            throw new IllegalArgumentException("Column " + column + " not found in file " + filename);
        }
        
        // Extract the values
        double[] doubles = new double[data.size() - 1]; // Exclude header
        for (int i = 1; i < data.size(); i++) {
            doubles[i - 1] = Double.parseDouble(data.get(i)[columnIndex]);
        }
        
        return round4(new Gamma(doubles, gamma, normalise).getSegmentPenalty());
    }
    
    /**
     * Calculates the gamma value for a file with default parameters.
     *
     * @param filename The name of the file
     * @return The gamma value
     */
    private double gammaForFile(String filename) {
        return gammaForFile(filename, 50.0, true, "ratio");
    }
}