package com.hartwig.hmftools.amber.contamination;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.amber.ImmutableBaseDepthData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

public class RealDataTest
{
    private File dataDir = new File("/Users/timlavers/work/batches/2026/1/30/3/data");
    private ChrArmLocator mChrArmLocator = ChrArmLocator.defaultLocator(RefGenomeVersion.V37);

    @Test
    public void heterozygosity()
    {
        List<Double> vafs = List.of(0.5);
        runSample("Sample_15838972", vafs);
        runSample("Sample_15395296", vafs);
        runSample("Sample_15418771", vafs);
        runSample("Sample_14822014", vafs);
        runSample("Sample_14916792", vafs);
    }

    @Test
    public void contamination()
    {
        // Estimated contamination 0.2850, but peaks seem to be at 0.16, 0.35, 0.5, 0.65, 0.84
        List<Double> vafs = new ArrayList<>();
        double vaf = 0.03;
        for(int i = 0; i < 92; i++)
        {
            vafs.add(vaf);
            vaf += 0.01;
        }
        runSample("Sample_15100025", vafs);
    }

    @Test
    public void sdTest()
    {
        //        List<Double> vafs = List.of(0.07, 0.14, 0.14, 0.5, 0.86, 0.93);
        List<Double> vafs = List.of(0.001, 0.002, 0.003);
        //        List<Integer> ns = List.of(10, 20, 50, 100, 200, 1000);
        List<Integer> ns = List.of(10, 100, 1000, 10000, 100000);
        for(double vaf : vafs)
        {
            System.out.println("VAF:" + vaf);
            for(int n : ns)
            {
                double sd = Math.sqrt(vaf * (1 - vaf) * n);
                double mean = vaf * n;
                int lower = (int) Math.floor(mean - 2 * sd);
                int upper = (int) Math.ceil(mean + 2 * sd);
                System.out.println(lower + " " + upper);
            }
        }
    }

    @Test
    public void exhaustiveSearchForPeaksTest()
    {
        List<Double> vafs = new ArrayList<>();
        double vaf = 0.01;
        for(int i = 0; i < 50; i++)
        {
            vafs.add(vaf);
            vaf += 0.01;
        }
        for(File dataFile : dataDir.listFiles((dir, name) -> name.endsWith("gz")))
        {
            List<VafResult> results = runSample(dataFile, vafs);
            List<VafResult> peaks = findPeaks(results);
            if(peaks.size() > 1)
            {
                System.out.println(dataFile.getName() + " " + peaks.size());
                System.out.println(peaks);
            }
        }

        //        List<VafResult> factors = runSample("Sample_15838972", vafs);
        //        List<VafResult> peaks = findPeaks(factors);
        //        System.out.println(peaks);
    }

    private List<VafResult> findPeaks(List<VafResult> factors)
    {
        List<VafResult> results = new ArrayList<>();
        VafResult previous = null;
        boolean descending = false;
        for(VafResult vafResult : factors)
        {
            VafConsistencyCheckResult factor = vafResult.check();
            if(factor.gini() > 0.3)
            {
                continue;
            }
            if(previous != null)
            {
                if(factor.totalWeightInBand() > previous.check().totalWeightInBand())
                {
                    descending = false;
                }
                if(factor.totalWeightInBand() == previous.check().totalWeightInBand())
                {
                    descending = false;
                }
                if(factor.totalWeightInBand() < previous.check().totalWeightInBand())
                {
                    if(!descending)
                    {
                        results.add(previous);
                    }
                    descending = true;
                }
            }
            previous = vafResult;
        }
        return results;
    }

    @Test
    public void contamination2()
    {
        // Estimated contamination 0.1, but peaks seem to be at 0.16, 0.35, 0.5, 0.65, 0.84
        List<Double> vafs = new ArrayList<>();
        double vaf = 0.001;
        for(int i = 0; i < 100; i++)
        {
            vafs.add(vaf);
            vaf += 0.001;
        }
        runSample("Sample_15837704", vafs);
        //        runSample("Sample_15100025", vafs);
        //        runSample("Sample_15100022", vafs);
    }

    @Test
    public void hcc1395()
    {
        File dataFile = new File("/Users/timlavers/work/batches/2026/2/4/1/HCC1395N.amber.raw.tsv.gz");
        List<TumorContamination> rawData = DelimFileReader.read(dataFile, new RowReader());
        List<TumorContamination> data = rawData
                .stream()
                .filter(t -> !t.chromosome().endsWith("X"))
                .filter(t -> !t.chromosome().endsWith("Y"))
                .filter(t -> t.Tumor.readDepth() > 20)
                .toList();

        int totalRefCount = data.stream().mapToInt(t -> t.Tumor.refSupport()).sum();
        int totalAltCount = data.stream().mapToInt(t -> t.Tumor.altSupport()).sum();

/*
        DescriptiveStatistics depthStats = new DescriptiveStatistics();
        data.forEach(t -> depthStats.addValue(t.Tumor.readDepth()));
        System.out.println(depthStats.getMin() + "-" + depthStats.getMax() + ", mean: " + depthStats.getMean() + ", median: "
                + depthStats.getPercentile(50) + ", sd: " + depthStats.getStandardDeviation());
        List<Double> vafs = new ArrayList<>();
        double vaf = 0.03;
        for(int i = 0; i < 92; i++)
        {
            vafs.add(vaf);
            vaf += 0.01;
        }
        for(double v : vafs)
        {
            var factor = PerArmVafConsistencyChecker.calculateConfirmationFactor(mChrArmLocator, v, data);
            printResult(v, factor);
        }
 */
    }

    @Test
    public void copyNumberEvents()
    {
        List<Double> vafs = List.of(0.07, 0.14, 0.14, 0.5, 0.86, 0.93);
        runSample("Sample_15020807", vafs);
    }

    private List<VafResult> runSample(final String sampleId, final List<Double> vafs)
    {
        System.out.println("\n" + sampleId);
        File sampleFile = new File(dataDir, sampleId + ".amber.raw.tsv.gz");
        return runSample(sampleFile, vafs);
    }

    private List<VafResult> runSample(File sampleFile, final List<Double> vafs)
    {
        List<TumorContamination> data = readSample(sampleFile)
                .stream()
                .filter(t -> !t.chromosome().endsWith("X"))
                .filter(t -> !t.chromosome().endsWith("Y"))
                .filter(t -> t.Tumor.readDepth() > 20)
                .toList();
        //        List<TumorContamination> data = readSample(sampleId);

        DescriptiveStatistics depthStats = new DescriptiveStatistics();
        data.forEach(t -> depthStats.addValue(t.Tumor.readDepth()));
        //        System.out.println(
        //                "Depth stats: " + depthStats.getMin() + "-" + depthStats.getMax() + ", mean: " + depthStats.getMean() + ", median: "
        //                        + depthStats.getPercentile(50) + ", sd: " + depthStats.getStandardDeviation());

        List<VafResult> results = new ArrayList<>();
        for(double vaf : vafs)
        {
            var factor = PerArmVafConsistencyChecker.calculateConfirmationFactor(mChrArmLocator, vaf, data);
            //            printResult(vaf, factor);
            results.add(new VafResult(vaf, factor));
        }
        return results;
    }

    void printResult(double vaf, VafConsistencyCheckResult factor)
    {
        System.out.printf("%.3f %.5f %d %d%n",
                vaf,
                factor.gini(),
                factor.totalWeightInBand(),
                factor.totalWeightAcrossAllVafValues());
    }

    private List<TumorContamination> readSample(File sampleFile)
    {
        return DelimFileReader.read(sampleFile, new RowReader());
    }
}

class RowReader implements Function<DelimFileReader.Row, TumorContamination>
{
    @Override
    public TumorContamination apply(DelimFileReader.Row row)
    {
        BaseDepthData tumorBDD = ImmutableBaseDepthData.builder()
                .ref(BaseDepthData.Base.valueOf(row.get(2)))
                .alt(BaseDepthData.Base.valueOf(row.get(3)))
                .readDepth(row.getInt(4))
                .indelCount(row.getInt(5))
                .refSupport(row.getInt(6))
                .altSupport(row.getInt(7))
                .build();
        return new TumorContamination(row.get(0), row.getInt(1), null, tumorBDD);
    }
}
