package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.linx.analyser.SvTestRoutines.createBnd;
import static com.hartwig.hmftools.linx.analyser.SvTestRoutines.createDel;
import static com.hartwig.hmftools.linx.analyser.SvTestRoutines.createDup;
import static com.hartwig.hmftools.linx.analyser.SvTestRoutines.createIns;
import static com.hartwig.hmftools.linx.analyser.SvTestRoutines.createInv;
import static com.hartwig.hmftools.linx.analyser.SvTestRoutines.createSgl;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.linx.analysis.SvUtilities;
import com.hartwig.hmftools.linx.stats.FisherExactTest;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

public class MiscTest
{

    @Test
    public void testConsistency()
    {
        final SvVarData del = createDel(1, "1", 100, 200);
        assertEquals(calcConsistency(del), 0);

        final SvVarData ins = createIns(1, "1", 100, 200);
        assertEquals(calcConsistency(ins), 0);

        final SvVarData dup = createDup(1, "1", 100, 200);
        assertEquals(calcConsistency(dup), 0);

        final SvVarData inv = createInv(1, "1", 100, 200, 1);
        assertEquals(calcConsistency(inv), 2);

        final SvVarData bnd = createBnd(1, "1", 100, 1, "2", 100, -1);
        assertEquals(calcConsistency(bnd), 0);

        final SvVarData sgl = createSgl(1, "1", 100, 1, false);
        assertEquals(calcConsistency(sgl), 1);
    }

    @Test
    public void testMiscMethods()
    {
        assertTrue(makeChrArmStr("1", "P").equals("1_P"));

        String test = "something";
        test = appendStr(test, "else", ';');

        assertEquals("something;else", test);
    }

    @Test
    public void testBreakendLists()
    {
        LinxTester tester = new LinxTester();

        // 4 breakends, 2 on each arm
        long centromerePos = SvUtilities.getChromosomalArmLength("1", CHROMOSOME_ARM_P);
        long qArmPos = centromerePos + 10000000;
        final SvVarData var1 = createDel(tester.nextVarId(), "1", 100,200);
        final SvVarData var2 = createDel(tester.nextVarId(), "1", 300,400);
        final SvVarData var3 = createDel(tester.nextVarId(), "1", qArmPos + 1000, qArmPos + 2000);
        final SvVarData var4 = createDel(tester.nextVarId(), "1", qArmPos + 10000,qArmPos + 20000);

        // add them out of order which will require partial chain reconciliation
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();

        List<SvBreakend> breakendList = tester.Analyser.getState().getChrBreakendMap().get("1");

        assertEquals(3, findCentromereBreakendIndex(breakendList, CHROMOSOME_ARM_P));
        assertEquals(4, findCentromereBreakendIndex(breakendList, CHROMOSOME_ARM_Q));

        // try again with breakends only in 1 list
        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();

        breakendList = tester.Analyser.getState().getChrBreakendMap().get("1");

        assertEquals(3, findCentromereBreakendIndex(breakendList, CHROMOSOME_ARM_P));
        assertEquals(-1, findCentromereBreakendIndex(breakendList, CHROMOSOME_ARM_Q));

        tester.clearClustersAndSVs();
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();

        breakendList = tester.Analyser.getState().getChrBreakendMap().get("1");

        assertEquals(-1, findCentromereBreakendIndex(breakendList, CHROMOSOME_ARM_P));
        assertEquals(0, findCentromereBreakendIndex(breakendList, CHROMOSOME_ARM_Q));
    }


    @Test
    public void testStatsRoutines()
    {
        FisherExactTest fetCalc = new FisherExactTest();
        fetCalc.initialise(1000);

        int withAwithB = 11;
        int withANoB = 27;
        int noAWithB = 2;
        int noAnoB = 170;
        double expectedCount = 5;

        double fisherProb = fetCalc.getLeftTailedP(withAwithB, noAWithB, withANoB, noAnoB);
        fisherProb = fetCalc.getRightTailedP(withAwithB, noAWithB, withANoB, noAnoB);

    }

}
