package com.hartwig.hmftools.amber.purity;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class CandidatePeakEvaluationTest
{
    private static final double LEVEL = 0.1;
    private CandidatePeak Peak;
    private int currentEvidencePosition;

    @Before
    public void setup()
    {
        Peak = new CandidatePeak(LEVEL);
        currentEvidencePosition = 1000;
    }

    @Test
    public void needToCaptureMinimumNumberOfDeepPointsTest()
    {
        // 15 points, but 5 are not deep enough.
        List<PositionEvidence> deepEnoughPoints = createHomPeakEvidencePointsWithDepth(10, 1000);
        List<PositionEvidence> shallowPoints = createHomPeakEvidencePointsWithDepth(5, 91); // See CandidatePeakTest
        List<PositionEvidence> allPoints = new ArrayList<>(deepEnoughPoints);
        allPoints.addAll(shallowPoints);
        CandidatePeakEvaluationResult result = calculateResult(allPoints);
        assertEquals(0.0, result.score(), 0.0001);

        // 15 deep enough points
        List<PositionEvidence> newPoints = createHomPeakEvidencePointsWithDepth(5, 920);
        List<PositionEvidence> enoughDeepPoints = new ArrayList<>(deepEnoughPoints);
        enoughDeepPoints.addAll(newPoints);
        result = calculateResult(enoughDeepPoints);
        Assert.assertTrue(result.score() > 0.0);
    }

    @Test
    public void needToCaptureAMinimumNumberOfPointsTest()
    {
        List<PositionEvidence> capturedPoints = createHomPeakEvidencePointsWithDepth(14, 1000);
        List<PositionEvidence> missedPoints = createUncapturedPointsWithDepth(14, 1000);
        List<PositionEvidence> allPoints = new ArrayList<>(capturedPoints);
        allPoints.addAll(missedPoints);
        CandidatePeakEvaluationResult result = calculateResult(allPoints);
        assertEquals(0.0, result.score(), 0.0001);

        List<PositionEvidence> moreCapturedPoints = createHomPeakEvidencePointsWithDepth(1, 1000);
        allPoints = new ArrayList<>(capturedPoints);
        allPoints.addAll(moreCapturedPoints);
        result = calculateResult(allPoints);
        Assert.assertTrue(result.score() > 0.0);
    }

    @Test
    public void capturedPointsCanBeInHetPeakTest()
    {
        List<PositionEvidence> points = createHetPeakEvidencePointsWithDepth(15, 1000);
        CandidatePeakEvaluation evaluation = new CandidatePeakEvaluation(Peak, points);
        CandidatePeakEvaluationResult result = calculateResult(points);
        Assert.assertTrue(result.score() > 0.0);
    }

    @Test
    public void scoreIsRatioOfCapturedPointsTest()
    {
        List<PositionEvidence> homPoints = createHomPeakEvidencePointsWithDepth(15, 1000);
        List<PositionEvidence> hetPoints = createHetPeakEvidencePointsWithDepth(16, 1000);
        List<PositionEvidence> uncapturedPoints = createUncapturedPointsWithDepth(69, 1000);
        List<PositionEvidence> allPoints = new ArrayList<>(homPoints);
        allPoints.addAll(hetPoints);
        allPoints.addAll(uncapturedPoints);
        CandidatePeakEvaluationResult result = calculateResult(allPoints);
        assertEquals(0.31, result.score(), 0.0001);
    }

    private CandidatePeakEvaluationResult calculateResult(List<PositionEvidence> evidencePoints)
    {
        CandidatePeakEvaluation evaluation = new CandidatePeakEvaluation(Peak, evidencePoints);
        evaluation.call();
        return evaluation.result();
    }

    private List<PositionEvidence> createHetPeakEvidencePointsWithDepth(int numberOfPoints, int depth)
    {
        int altDepth = (int) (LEVEL / 2.0 * depth);
        int refDepth = depth - altDepth;
        return createPoints(numberOfPoints, depth, refDepth, altDepth);
    }

    private List<PositionEvidence> createHomPeakEvidencePointsWithDepth(int numberOfPoints, int depth)
    {
        int altDepth = (int) (LEVEL * depth);
        int refDepth = depth - altDepth;
        return createPoints(numberOfPoints, depth, refDepth, altDepth);
    }

    private List<PositionEvidence> createUncapturedPointsWithDepth(int numberOfPoints, int depth)
    {
        return createPoints(numberOfPoints, depth, depth, 0);
    }

    @NotNull
    private List<PositionEvidence> createPoints(final int numberOfPoints, final int depth, final int refDepth, final int altDepth)
    {
        List<PositionEvidence> evidencePoints = new ArrayList<>(numberOfPoints);
        for(int i = 0; i < numberOfPoints; i++)
        {
            PositionEvidence evidencePoint = new PositionEvidence("1", currentEvidencePosition, "A", "T");
            evidencePoint.ReadDepth = depth;
            evidencePoint.RefSupport = refDepth;
            evidencePoint.AltSupport = altDepth;
            evidencePoints.add(evidencePoint);
            currentEvidencePosition += 100;
        }
        return evidencePoints;
    }
}
