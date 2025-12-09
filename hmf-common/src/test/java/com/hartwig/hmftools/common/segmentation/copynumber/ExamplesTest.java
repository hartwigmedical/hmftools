package com.hartwig.hmftools.common.segmentation.copynumber;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class ExamplesTest extends SegmentationTestBase
{
    @Test
    public void amber()
    {
        PiecewiseConstantFit pcf = pcfForFile("amber1.tsv", 100.0, true, "tumorModifiedBAF");
        PiecewiseConstantFit expected = new PiecewiseConstantFit(new int[] { 5, 5, 5 }, new int[] { 0, 5, 10 }, d(0.515, 0.962, 0.511));
        assertEquals(expected, pcf);
    }

    @Test
    public void amber2()
    {
        double[] rounded = readDoubles("amber2.tsv", "tumorModifiedBAF");
        final double penalty = 0.036533517458912736;
        PenaltyCalculator penaltyCalculator = new FixedPenalty(penalty);
        Segmenter calculation = new Segmenter(rounded, penaltyCalculator);
        // calculation.leastCostSegmentation.cost(0.036533517458912736) = 0.225
        //        Segmentation fromR = new Segmentation()
        List<double[]> rIntervals = splitIntoListOfSubArraysOfGivenLengths(List.of(1, 1, 108, 1, 19), rounded);
        Segmentation rResult = new Segmentation(rIntervals);
        double rCost = rResult.cost(penalty);
        Assert.assertTrue(calculation.leastCostSegmentation.cost(penalty) < rCost);
    }

    private static List<double[]> splitIntoListOfSubArraysOfGivenLengths(List<Integer> lengths, double[] doubles)
    {
        List<double[]> result = new ArrayList<>();
        int cursor = 0;
        for(int length : lengths)
        {
            result.add(Arrays.copyOfRange(doubles, cursor, cursor + length));
            cursor += length;
        }
        return result;
    }

    @Test
    public void minima()
    {
        assertEquals(
                new PiecewiseConstantFit(new int[] { 16 }, new int[] { 0 }, d(0.883)),
                pcfForFile("ratios1.tsv")
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 36 }, new int[] { 0 }, d(0.902)),
                pcfForFile("ratios2.tsv")
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 75 }, new int[] { 0 }, d(0.899)),
                pcfForFile("ratios3.tsv")
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 109 }, new int[] { 0 }, d(0.889)),
                pcfForFile("ratios4.tsv")
        );
    }

    @Test
    public void panel9pSample1()
    {
        String filename = "ratios_tii_9p.tsv"; // COREDB010001TII
        // Note that the R code gives the wrong mean for the first segment. It gives -0.156017.
        assertEquals(
                new PiecewiseConstantFit(new int[] { 86, 4, 19 }, new int[] { 0, 86, 90 }, d(-0.153, -0.764, -0.193)),
                pcfForFile(filename, 35.0, true, "logRatio")
        );
        // Note that the R code incorrectly has a length-1 segment at the beginning.
        assertEquals(
                new PiecewiseConstantFit(
                        new int[] { 7, 2, 77, 4, 3, 16 },
                        new int[] { 0, 7, 9, 86, 90, 93 },
                        d(-0.159, -0.594, -0.141, -0.764, 0.123, -0.252)
                ),
                pcfForFile(filename, 10.0, true, "logRatio")
        );
    }

    @Test
    public void panel9pSample2()
    {
        String filename = "ratios_27t_9p.tsv"; // COREDB011427T
        // Assertions made by comparing with the output from the R code and excel (as the R code
        // always adds a spurious length 1 segment).
        assertEquals(
                new PiecewiseConstantFit(new int[] { 109 }, new int[] { 0 }, d(-0.227)),
                pcfForFile(filename, 50.0, true, "logRatio")
        );
        assertEquals(
                new PiecewiseConstantFit(new int[] { 86, 7, 16 }, new int[] { 0, 86, 93 }, d(-0.195, -0.604, -0.236)),
                pcfForFile(filename, 18.0, true, "logRatio")
        );
    }

    @Test
    public void chr1p()
    {
        // This takes about a minute to run as there are 108000 points in the file.
        // This is about as large a file as we can expect to process.
        // Processing of this input failed with an early version of the code.
        String filename = "chr1p.ratios.tsv";
        /*
        R output:
        $Lengde
  [1]    1    1    5    4    3   54    1  250    3  516    1  700  351
 [14]    1  609    1  700    1    3    6  785  650    1  157    1    1
 [27] 1768    1  707    2 2407    7 2259    3  572   11    3    1   14
 [40]    3    8    5   19 1044    1 6726    1  680    1 2792 1411    1
 [53]  804    1  263 2665   10 1452 1642    1   46  712    1  348 3842
 [66]    1 1488    1 1762 4127 4064    2 1188 1310    1 1455 1676    1
 [79]  548    1  575  650  786 1767 1254    1  594    1 2002    1  732
 [92] 1180   46 1550    1 1758    2 1375    1 1274 2080    3  416    1
[105] 2405    4 1223   26 2629 2784    1  355    1  922    1  992    3
[118] 1246 5007  715    1 5671    1  583    4  570    2  570    1    1
[131]   23    4  713  349    9 1175    1   43    6    9    1 1136    1
[144] 2767 3367    2    1    1

Our's agrees with this excep that the first two length-1 segments
given by R are incorrect (have a look at the data).
         */

        PiecewiseConstantFit calculated = pcfForFile(filename, 50.0, true, "logRatio");
        assertEquals(146, calculated.lengths().length);
        assertEquals(7, calculated.lengths()[0]);
        assertEquals(4, calculated.lengths()[1]);
        assertEquals(3, calculated.lengths()[2]);
        assertEquals(54, calculated.lengths()[3]);

        assertEquals(2767, calculated.lengths()[141]);
        assertEquals(2, calculated.lengths()[143]);
        assertEquals(1, calculated.lengths()[145]);
    }

    @Test
    public void chr2q()
    {
        // Chr2.q is the longest chromosome arm in the human genome.
        // There are about 135,000 points in the input file.
        // This data is from COLO829T.
        String filename = "chr2q.ratios.tsv";
        /*
        R output:
        segmentPenalty: 0.7133056 (we get 0.71212)

Note the usual spurious segment at the start...
$Lengde
 [1]    1 5975    1 1032    1 1559 1163    5 2016    1 3103 3638 3227    1
[15] 4620    2  978    2 3534    2 2112 7969 3058   14 9943 2263    1 9790
[29]    2 2557    2 6674   11 3625 5057 3734 3900    1 5275 1130 7354    2
[43] 3973 1736    2 2015    3 2618 2394 2029    1  864   18 2378  108    3
[57] 3169 6964    1 1225

$sta
 [1]      1      2   5977   5978   7010   7011   8570   9733   9738  11754
[11]  11755  14858  18496  21723  21724  26344  26346  27324  27326  30860
[21]  30862  32974  40943  44001  44015  53958  56221  56222  66012  66014
[31]  68571  68573  75247  75258  78883  83940  87674  91574  91575  96850
[41]  97980 105334 105336 109309 111045 111047 113062 113065 115683 118077
[51] 120106 120107 120971 120989 123367 123475 123478 126647 133611 133612         */

        PiecewiseConstantFit calculated = pcfForFile(filename, 100.0, true, "logRatio");
        assertEquals(59, calculated.lengths().length);
        assertEquals(5976, calculated.lengths()[0]);
        assertEquals(1, calculated.lengths()[1]);
        assertEquals(1032, calculated.lengths()[2]);
        assertEquals(1, calculated.lengths()[3]);
        assertEquals(1559, calculated.lengths()[4]);

        assertEquals(59, calculated.startPositions().length);
        assertEquals(0, calculated.startPositions()[0]);
        assertEquals(5976, calculated.startPositions()[1]);
        assertEquals(5977, calculated.startPositions()[2]);
        assertEquals(7009, calculated.startPositions()[3]);
        assertEquals(7010, calculated.startPositions()[4]);
    }

    @Test
    public void chr10p()
    {
        String filename = "chr10p.ratios.tsv";
        double[] doubles = readDoubles(filename, "logRatio");
        /*
        The segmentation we calculate is:
[354, 517, 1897, 313, 828, 411, 2, 1022, 3, 575, 2115, 1875, 5092, 671, 1275, 1461, 1785, 2, 2099, 1073, ...and 48 more (set the 'kotest.assertions.collection.print.size' JVM property to see more / less items)]>
This is close to but not the same as calculated by R.
The cost of the two solutions is the same though.
         */
        Segmenter calculation = new Segmenter(doubles, 50.0, true);
        Segmentation rResult = buildSegmentationForROutput(doubles);
        double calculatedCost = rResult.cost(calculation.pcf().means()[0]); // Using the first mean as a proxy for segmentPenalty
        double rCost = rResult.cost(calculation.pcf().means()[0]);
        assertEquals(calculatedCost, rCost, 1e-10);

        assertEquals(
                new PiecewiseConstantFit(
                        new int[] { 35543, 89, 51 },
                        new int[] { 0, 35543, 35632 },
                        d(-0.007, 1.394, 2.459)
                ),
                pcfForFile(filename, 2000.0, true, "logRatio")
        );
    }

    @Test
    public void chr10q()
    {
        String filename = "chr10q.ratios.tsv";
        assertEquals(
                new PiecewiseConstantFit(
                        new int[] { 24, 104, 5, 39, 2849, 174, 72458, 15, 6684 },
                        new int[] { 0, 24, 128, 133, 172, 3021, 3195, 75653, 75668 },
                        d(5.135, 0.264, 6.383, 1.076, -0.021, 0.961, -0.016, 2.88, -0.016)
                ),
                pcfForFile(filename, 5000.0, true, "logRatio")
        );
    }

    @Test
    public void chr13()
    {
        assertEquals(
                new PiecewiseConstantFit(new int[] { 108, 25, 48 }, new int[] { 0, 108, 133 }, d(0.546, 0.185, 0.519)),
                pcfForFile("ratios_chr13_1.tsv", 50.0, true, "ratio")
        );
    }

    private Segmentation buildSegmentationForROutput(double[] ratios)
    {
        List<Integer> intervalLengths = extractIntervalLengths();
        int sum = 0;
        for(int length : intervalLengths)
        {
            sum += length;
        }
        assertEquals(ratios.length, sum);

        List<double[]> intervals = new ArrayList<>();
        int cursor = 0;
        for(int intervalLength : intervalLengths)
        {
            int stop = cursor + intervalLength;
            double[] interval = new double[intervalLength];
            System.arraycopy(ratios, cursor, interval, 0, intervalLength);
            intervals.add(interval);
            cursor = stop;
        }
        return new Segmentation(intervals);
    }

    private List<Integer> extractIntervalLengths()
    {
        List<Integer> integers = new ArrayList<>();
        try(BufferedReader reader = new BufferedReader(new FileReader(loadFileFromResources("r_output_1.txt"))))
        {
            String line;
            while((line = reader.readLine()) != null)
            {
                String[] parts = line.split(",");
                for(String part : parts)
                {
                    integers.add(Integer.parseInt(part.trim()));
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("Error reading file: " + "r_output_1.txt", e);
        }
        return integers;
    }

    private PiecewiseConstantFit pcfForFile(String filename, double gamma, boolean normalise, String column)
    {
        double[] rounded = readDoubles(filename, column);
        return new Segmenter(rounded, gamma, normalise).pcf();
    }

    private PiecewiseConstantFit pcfForFile(String filename)
    {
        return pcfForFile(filename, 50.0, false, "ratio");
    }

    private double[] readDoubles(String filename, String column)
    {
        List<String[]> data = readResourcesFile(filename);

        // Find the column index
        int columnIndex = -1;
        if(!data.isEmpty())
        {
            String[] header = data.get(0);
            for(int i = 0; i < header.length; i++)
            {
                if(header[i].equals(column))
                {
                    columnIndex = i;
                    break;
                }
            }
        }

        if(columnIndex == -1)
        {
            throw new IllegalArgumentException("Column " + column + " not found in file " + filename);
        }

        // Extract the values
        double[] doubles = new double[data.size() - 1]; // Exclude header
        for(int i = 1; i < data.size(); i++)
        {
            doubles[i - 1] = Double.parseDouble(data.get(i)[columnIndex]);
        }

        // Round the values
        double[] rounded = new double[doubles.length];
        for(int i = 0; i < doubles.length; i++)
        {
            rounded[i] = round10(doubles[i]);
        }

        System.out.println("Number of data points: " + rounded.length);
        return rounded;
    }
}
