package com.hartwig.hmftools.linx.misc;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDel;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createDup;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createIns;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createInv;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSgl;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.linx.analysis.SvUtilities;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.utils.LinxTester;

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

        final SvVarData sgl = createSgl(1, "1", 100, 1);
        assertEquals(calcConsistency(sgl), 1);
    }

    @Test
    public void testBreakendLists()
    {
        LinxTester tester = new LinxTester();

        // 4 breakends, 2 on each arm
        int centromerePos = SvUtilities.getChromosomalArmLength("1", P_ARM);
        int qArmPos = centromerePos + 10000000;
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

        assertEquals(3, findCentromereBreakendIndex(breakendList, P_ARM));
        assertEquals(4, findCentromereBreakendIndex(breakendList, Q_ARM));

        // try again with breakends only in 1 list
        tester.clearClustersAndSVs();
        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);

        tester.preClusteringInit();

        breakendList = tester.Analyser.getState().getChrBreakendMap().get("1");

        assertEquals(3, findCentromereBreakendIndex(breakendList, P_ARM));
        assertEquals(-1, findCentromereBreakendIndex(breakendList, Q_ARM));

        tester.clearClustersAndSVs();
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();

        breakendList = tester.Analyser.getState().getChrBreakendMap().get("1");

        assertEquals(-1, findCentromereBreakendIndex(breakendList, P_ARM));
        assertEquals(0, findCentromereBreakendIndex(breakendList, Q_ARM));
    }
}
