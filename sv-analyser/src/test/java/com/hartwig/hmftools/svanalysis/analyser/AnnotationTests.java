package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createBnd;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createInv;
import static com.hartwig.hmftools.svanalysis.analyser.SvTestHelper.createSgl;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.junit.Test;

public class AnnotationTests
{
    @Test
    public void testUnderClusteredFoldbackAnnotations()
    {
        SvTestHelper tester = new SvTestHelper();
        tester.logVerbose(true);

        double[] chrCopyNumbers = {2.0, 2.0, 2.0};
        tester.ClusteringMethods.getChrCopyNumberMap().put("1", chrCopyNumbers);
        tester.ClusteringMethods.getChrCopyNumberMap().put("2", chrCopyNumbers);

        final SvVarData var1 = createInv("0", "1", 101000, 104000, -1);

        // straddling BNDs
        final SvVarData var2 = createBnd("1", "1", 1000, 1, "2", 1000, -1);
        final SvVarData var3 = createBnd("2", "1", 200000, 1, "2", 1100, 1);

        // single other cluster
        final SvVarData var4 = createSgl("3", "1", 150000, 1, false);

        tester.AllVariants.add(var1);
        tester.AllVariants.add(var2);
        tester.AllVariants.add(var3);
        tester.AllVariants.add(var4);

        tester.preClusteringInit();

        tester.Analyser.clusterAndAnalyse();

        assertEquals(tester.Analyser.getClusters().size(), 3);

    }
}
