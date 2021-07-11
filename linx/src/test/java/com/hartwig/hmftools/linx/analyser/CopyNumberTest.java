package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createTestSv;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.populateChromosomeBreakendMap;
import static com.hartwig.hmftools.linx.cn.CnJcnCalcs.calcAdjustedJcnValues;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.List;

import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class CopyNumberTest
{
    @Test
    public void testPloidyRecalcs()
    {
        // example
        double cnChgStart = 3.1;
        double cnChgEnd = 3.6;
        int tumorReadCount = 58;
        double maxCNStart = 6.01;
        double maxCNEnd = 6.66;
        double ploidy = 2.14;
        int[] startDepthData = {3, 3};
        int[] endDepthData = {1, 1};

        JcnCalcData calcResults =  calcAdjustedJcnValues(cnChgStart, cnChgEnd, tumorReadCount, ploidy,
                maxCNStart, maxCNEnd, startDepthData, endDepthData);

        assertEquals(2.93, calcResults.JcnEstimate, 0.01);
        assertEquals(0.83, calcResults.JcnUncertainty, 0.01);

        cnChgStart = 2;
        cnChgEnd = 1.3;
        tumorReadCount = 69;
        maxCNStart = 7;
        maxCNEnd = 5;
        ploidy = 1.43;
        startDepthData[0] = 0;;
        endDepthData[0] = 0;

        calcResults =  calcAdjustedJcnValues(cnChgStart, cnChgEnd, tumorReadCount, ploidy,
                maxCNStart, maxCNEnd, startDepthData, endDepthData);

        assertEquals(1.48, calcResults.JcnEstimate, 0.01);
        assertEquals(1.69, calcResults.JcnUncertainty, 0.01);
    }

    @Test
    public void testClusterMinMaxCopyNumber()
    {
        LinxTester tester = new LinxTester();
        // tester.logVerbose(true);

        SvVarData var1 = createDel(tester.nextVarId(), "1", 100, 200);
        SvVarData var2 = createDel(tester.nextVarId(), "1", 300, 400);
        SvVarData var3 = createDel(tester.nextVarId(), "1", 500, 600);

        var1.setJcnRecalcData(0.9, 1.6);
        var2.setJcnRecalcData(0.1, 1.2);
        var3.setJcnRecalcData(0.99, 2.01);

        SvCluster cluster1 = new SvCluster(0);
        cluster1.addVariant(var1);
        cluster1.addVariant(var2);
        cluster1.addVariant(var3);

        assertFalse(cluster1.hasVariedJcn());
        assertEquals(1.0, cluster1.getJcnRange()[SE_START], 0.001);
        assertEquals(1.0, cluster1.getJcnRange()[SE_END], 0.001);

        // now check alignment to a common but non-integer ploidy
        SvVarData var4 = createDel(tester.nextVarId(), "1", 700, 800);

        // all SVs have 1.1-1.2 in common
        var4.setJcnRecalcData(1.1, 1.5);
        cluster1.addVariant(var4);

        assertFalse(cluster1.hasVariedJcn());
        assertEquals(1.0, cluster1.getJcnRange()[SE_START], 0.001);
        assertEquals(1.0, cluster1.getJcnRange()[SE_END], 0.001);
    }

    private void buildCopyNumberData(LinxTester tester)
    {
        tester.Analyser.getState().reset();
        populateChromosomeBreakendMap(tester.AllVariants, tester.Analyser.getState());
        tester.populateCopyNumberData(false);
    }

    @Test
    public void testCopyNumberReconstruction()
    {
        // test copy number data re-creation from SVs
        LinxTester tester = new LinxTester();
        // tester.logVerbose(true);

        final String chromosome = "1";

        // simple DEL, nothing else
        SvVarData var1 = createTestSv(1, chromosome, chromosome, 1000,2000, 1, -1, DEL,  1);

        tester.AllVariants.add(var1);
        buildCopyNumberData(tester);

        List<SvCNData> cnDataList = tester.CnDataLoader.getChrCnDataMap().get(chromosome);
        assertEquals(4, cnDataList.size());

        double delta = 0.01;

        assertEquals(2, cnDataList.get(0).CopyNumber, delta);
        assertEquals(0.5, cnDataList.get(0).ActualBaf, delta);

        assertEquals(1, cnDataList.get(1).CopyNumber, delta);
        assertEquals(1, cnDataList.get(1).ActualBaf, delta);

        assertEquals(2, cnDataList.get(2).CopyNumber, delta);

        // simple DUP with ploidy 2
        tester.clearClustersAndSVs();
        var1 = createTestSv(1, chromosome, chromosome, 1000,2000, -1, 1, DUP,  2);

        tester.AllVariants.add(var1);
        buildCopyNumberData(tester);

        cnDataList = tester.CnDataLoader.getChrCnDataMap().get(chromosome);
        assertEquals(4, cnDataList.size());

        assertEquals(3, cnDataList.get(0).CopyNumber, delta);
        assertEquals(0.67, cnDataList.get(0).ActualBaf, delta);

        assertEquals(5, cnDataList.get(1).CopyNumber, delta);
        assertEquals(0.8, cnDataList.get(1).ActualBaf, delta);

        assertEquals(3, cnDataList.get(2).CopyNumber, delta);

        // foldbacks facing away from P-arm telomere
        tester.clearClustersAndSVs();
        var1 = createTestSv(1, chromosome, chromosome, 1000,2000, -1, -1, INV,  2);
        SvVarData var2 = createTestSv(2, chromosome, chromosome, 5000,6000, 1, 1, INV,  1);
        SvVarData var3 = createTestSv(3, chromosome, "2", 8000,100, 1, 1, BND,  1);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        buildCopyNumberData(tester);

        cnDataList = tester.CnDataLoader.getChrCnDataMap().get(chromosome);
        assertEquals(7, cnDataList.size());

        assertEquals(1, cnDataList.get(0).CopyNumber, delta);
        assertEquals(1.00, cnDataList.get(0).ActualBaf, delta);

        assertEquals(3, cnDataList.get(1).CopyNumber, delta);
        assertEquals(0.67, cnDataList.get(1).ActualBaf, delta);

        assertEquals(5, cnDataList.get(2).CopyNumber, delta);
        assertEquals(0.8, cnDataList.get(2).ActualBaf, delta);

        assertEquals(4, cnDataList.get(3).CopyNumber, delta);
        assertEquals(0.75, cnDataList.get(3).ActualBaf, delta);

    }



}
