package com.hartwig.hmftools.common.blastn;

import com.hartwig.hmftools.common.genome.region.Strand;

import org.jetbrains.annotations.NotNull;

/**
 * Represents a BLASTn match with various alignment details.
 */
public class BlastnMatch
{
    // Fields
    private final int querySeqLen;
    @NotNull private final String subjectTitle;
    private final double percentageIdent;
    private final double queryCoverage;
    private final int alignmentLength;
    private final int numMismatch;
    private final int numGapOpenings;
    private int queryAlignStart;
    private int queryAlignEnd;
    private final int subjectAlignStart;
    private final int subjectAlignEnd;
    @NotNull private final Strand subjectFrame;
    private final double expectedValue;
    private final double bitScore;
    @NotNull private final String alignedPartOfQuerySeq;
    @NotNull private final String alignedPartOfSubjectSeq;

    public BlastnMatch(int querySeqLen, @NotNull String subjectTitle, double percentageIdent, double queryCoverage, int alignmentLength,
            int numMismatch, int numGapOpenings, int queryAlignStart, int queryAlignEnd, int subjectAlignStart, int subjectAlignEnd,
            @NotNull Strand subjectFrame, double expectedValue, double bitScore, @NotNull String alignedPartOfQuerySeq,
            @NotNull String alignedPartOfSubjectSeq)
    {
        this.querySeqLen = querySeqLen;
        this.subjectTitle = subjectTitle;
        this.percentageIdent = percentageIdent;
        this.queryCoverage = queryCoverage;
        this.alignmentLength = alignmentLength;
        this.numMismatch = numMismatch;
        this.numGapOpenings = numGapOpenings;
        this.queryAlignStart = queryAlignStart;
        this.queryAlignEnd = queryAlignEnd;
        this.subjectAlignStart = subjectAlignStart;
        this.subjectAlignEnd = subjectAlignEnd;
        this.subjectFrame = subjectFrame;
        this.expectedValue = expectedValue;
        this.bitScore = bitScore;
        this.alignedPartOfQuerySeq = alignedPartOfQuerySeq;
        this.alignedPartOfSubjectSeq = alignedPartOfSubjectSeq;
    }

    public int getQuerySeqLen()
    {
        return querySeqLen;
    }

    @NotNull
    public String getSubjectTitle()
    {
        return subjectTitle;
    }

    public double getPercentageIdent()
    {
        return percentageIdent;
    }

    public double getQueryCoverage()
    {
        return queryCoverage;
    }

    public int getAlignmentLength()
    {
        return alignmentLength;
    }

    public int getNumMismatch()
    {
        return numMismatch;
    }

    public int getNumGapOpenings()
    {
        return numGapOpenings;
    }

    public int getQueryAlignStart()
    {
        return queryAlignStart;
    }

    public void setQueryAlignStart(final int queryAlignStart)
    {
        this.queryAlignStart = queryAlignStart;
    }

    public int getQueryAlignEnd()
    {
        return queryAlignEnd;
    }

    public void setQueryAlignEnd(final int queryAlignEnd)
    {
        this.queryAlignEnd = queryAlignEnd;
    }

    public int getSubjectAlignStart()
    {
        return subjectAlignStart;
    }

    public int getSubjectAlignEnd()
    {
        return subjectAlignEnd;
    }

    @NotNull
    public Strand getSubjectFrame()
    {
        return subjectFrame;
    }

    public double getExpectedValue()
    {
        return expectedValue;
    }

    public double getBitScore()
    {
        return bitScore;
    }

    @NotNull
    public String getAlignedPartOfQuerySeq()
    {
        return alignedPartOfQuerySeq;
    }

    @NotNull
    public String getAlignedPartOfSubjectSeq()
    {
        return alignedPartOfSubjectSeq;
    }
}
