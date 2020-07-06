package com.hartwig.hmftools.linx.analyser;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createBnd;
import static com.hartwig.hmftools.linx.utils.SvTestUtils.createSv;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.KNOWN_LINE_ELEMENT;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.POLY_A_MOTIF;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.SUSPECTED_LINE_ELEMENT;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.junit.Test;

import com.hartwig.hmftools.linx.utils.LinxTester;

public class LineTest
{

    @Test
    public void testKnownLineMarking()
    {
        LinxTester tester = new LinxTester();

        // scenario 1:
        // 2 BNDs within 5KB and 1 have poly A/T and not forming a DB
        // if other breakends form a DB, don't mark

        SvVarData bnd1 = createBnd(tester.nextVarId(), "1", 100, -1, "2", 100, 1);
        SvVarData bnd2 = createBnd(tester.nextVarId(), "1", 200, 1, "2", 200, -1);

        bnd1.setLineElement(KNOWN_LINE_ELEMENT, true);
        bnd2.setLineElement(KNOWN_LINE_ELEMENT, true);

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);

        tester.addClusterAndSVs(cluster);

        LineElementAnnotator leAnnotator = new LineElementAnnotator();

        int proximity = 5000;

        leAnnotator.markLineCluster(cluster, proximity);

        assertTrue(cluster.hasLinkingLineElements());
    }

    @Test
    public void testSuspectLineMarking()
    {
        LinxTester tester = new LinxTester();

        LineElementAnnotator leAnnotator = new LineElementAnnotator();
        int proximity = 5000;

        // scenario 1:
        // 2 BNDs within 5KB and 1 have poly A/T and not forming a DB with other breakends forming a short DB
        SvVarData bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 150,  -1, 1, BND, POLY_A_MOTIF);
        SvVarData bnd2 = createSv(tester.nextVarId(), "1", "2", 200, 170,  1, -1, BND, "");

        SvCluster cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        // in this case because the DB isn't short on the originating arm, both ends are marked as suspect
        assertEquals(bnd1.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertEquals(bnd1.getLineElement(false), NO_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(false), NO_LINE_ELEMENT);
        assertTrue(cluster.hasLinkingLineElements());

        // now with BNDs going to different arms
        tester.clearClustersAndSVs();

        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, POLY_A_MOTIF);
        bnd2 = createSv(tester.nextVarId(), "1", "3", 200, 110,  1, -1, BND, "");

        cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertTrue(bnd1.getLineElement(true).equals(SUSPECTED_LINE_ELEMENT));
        assertTrue(bnd2.getLineElement(true).equals(SUSPECTED_LINE_ELEMENT));
        assertEquals(bnd1.getLineElement(false), NO_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(false), NO_LINE_ELEMENT);
        assertTrue(cluster.hasLinkingLineElements());

        // now a BND with a remote SGL in a short DB
        tester.clearClustersAndSVs();

        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, POLY_A_MOTIF);
        SvVarData sgl = createSv(tester.nextVarId(), "2", "", 90, -1,  -1, -1, SGL, "");

        cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(sgl);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertTrue(bnd1.getLineElement(true).equals(SUSPECTED_LINE_ELEMENT));
        assertTrue(sgl.getLineElement(true).equals(NO_LINE_ELEMENT));
        assertEquals(bnd1.getLineElement(false), NO_LINE_ELEMENT);
        assertTrue(cluster.hasLinkingLineElements());

        // now test BNDs in a DB on the LINE arm which will invalidate the line test
        tester.clearClustersAndSVs();

        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "1", "2", 110, 110,  -1, -1, BND, POLY_A_MOTIF);

        cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertTrue(bnd1.getLineElement(true).equals(NO_LINE_ELEMENT));
        assertTrue(bnd2.getLineElement(true).equals(NO_LINE_ELEMENT));
        assertEquals(bnd1.getLineElement(false), NO_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(false), NO_LINE_ELEMENT);
        assertTrue(!cluster.hasLinkingLineElements());

        // now test BNDs in a longer DB on the remote arm which will invalidate the line test
        tester.clearClustersAndSVs();

        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "1", "2", 200, 200,  1, -1, BND, POLY_A_MOTIF);

        cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertTrue(bnd1.getLineElement(true).equals(NO_LINE_ELEMENT));
        assertTrue(bnd2.getLineElement(true).equals(NO_LINE_ELEMENT));
        assertEquals(bnd1.getLineElement(false), NO_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(false), NO_LINE_ELEMENT);
        assertTrue(!cluster.hasLinkingLineElements());

        // some other variant having the poly A motif
        tester.clearClustersAndSVs();

        bnd1 = createSv(tester.nextVarId(), "1", "2", 100, 100,  -1, 1, BND, "");
        bnd2 = createSv(tester.nextVarId(), "1", "3", 150, 200,  -1, -1, BND, "");
        SvVarData del = createSv(tester.nextVarId(), "1", "1", 1000, 11000,  1, -1, DEL, POLY_A_MOTIF);

        cluster = new SvCluster(0);
        cluster.addVariant(bnd1);
        cluster.addVariant(bnd2);
        cluster.addVariant(del);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertEquals(bnd1.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertEquals(del.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertEquals(del.getLineElement(false), NO_LINE_ELEMENT);
        assertEquals(bnd1.getLineElement(false), NO_LINE_ELEMENT);
        assertEquals(bnd2.getLineElement(false), NO_LINE_ELEMENT);
        assertTrue(cluster.hasLinkingLineElements());

        // 2 proximate breakends with poly A or T - only mark the breakends as line if they have the same orientation
        tester.clearClustersAndSVs();

        SvVarData sgl1 = createSv(tester.nextVarId(), "1", "0", 1000, 0,  -1, 0, SGL, POLY_A_MOTIF);
        SvVarData sgl2 = createSv(tester.nextVarId(), "1", "0", 1020, 0,  1, 0, SGL, POLY_A_MOTIF);

        cluster = new SvCluster(0);
        cluster.addVariant(sgl1);
        cluster.addVariant(sgl2);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertEquals(sgl1.getLineElement(true), NO_LINE_ELEMENT);
        assertEquals(sgl2.getLineElement(true), NO_LINE_ELEMENT);
        assertTrue(cluster.hasLinkingLineElements());

        tester.clearClustersAndSVs();

        sgl1 = createSv(tester.nextVarId(), "1", "0", 1000, 0,  -1, 0, SGL, POLY_A_MOTIF);
        sgl2 = createSv(tester.nextVarId(), "1", "0", 1100, 0,  -1, 0, SGL, POLY_A_MOTIF);

        cluster = new SvCluster(0);
        cluster.addVariant(sgl1);
        cluster.addVariant(sgl2);

        tester.addClusterAndSVs(cluster);
        tester.preClusteringInit();

        leAnnotator.markLineCluster(cluster, proximity);

        assertEquals(sgl1.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertEquals(sgl2.getLineElement(true), SUSPECTED_LINE_ELEMENT);
        assertTrue(cluster.hasLinkingLineElements());
    }

}