package com.hartwig.hmftools.telo.breakend;

import static com.hartwig.hmftools.telo.TeloConfig.TE_LOGGER;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.samtools.SamRecordUtils;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.telo.TeloConstants;
import com.hartwig.hmftools.telo.TeloUtils;
import com.hartwig.hmftools.telo.breakend.TelomericBreakEnd;
import com.hartwig.hmftools.telo.util.TelomereMatcher;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

// go through the bam file and find any potential
// split read sites
public class TelomericSplitReadAnalyser
{
    private class CandidateBreakEnd
    {
        public TelomericBreakEnd.Type Type;
        public String Chromosome;
        public int FirstPosition;
        public int LastPosition;
        public int NumSupportingFragments = 0;
    }

    // we let some bases at the edge to be non telomeric
    private double mTelomereMatchThreshold = 0.9;
    private int mBoundaryZone = 1;
    private int mMinSoftClipLength = 10;

    // the distance which we deem the same breakpoint
    private int mFuzzyMatchDistance = 10;

    private int mMinSupportingFragments = 1;

    private List<TelomericBreakEnd> mPotentialBreakEnds = new ArrayList<>();
    private List<TelomericBreakEnd> mConsolidatedBreakEnds = new ArrayList<>();

    void setBoundaryZone(int bz)
    {
        mBoundaryZone = bz;
    }

    void setMinSoftClipLength(int l)
    {
        mMinSoftClipLength = l;
    }

    void setFuzzyMatchDistance(int d) { mFuzzyMatchDistance = d; }

    // we want to collect evidence to support position of a break end
    public void processRead(@NotNull final SAMRecord r)
    {
        if (r.getReferenceName().equals("17") && r.getAlignmentStart() <= 22208891 && r.getAlignmentEnd() >= 22208891)
        {
            TE_LOGGER.info("here");
        }

        TelomericBreakEnd tbe = findTelomericBreakEnd(r);

        if (tbe != null)
        {
            mPotentialBreakEnds.add(tbe);
            //TE_LOGGER.trace("record: {} breakend: {} cigar: {} readBases: {}", r, tbe, r.getCigarString(), r.getReadString());
        }
    }

    // we want to collect all these potential break ends and clean them up into distinct
    // locations, and then find supporting or nonsupporting evidence
    public void consolidatePotentialBreakEnds()
    {
        List<CandidateBreakEnd> consolidatedBreakEnds = new ArrayList<>();

        // first sort all the breakends
        // we want to sort it such that position is sorted last, and then we can use that
        // to compare if positions are close by, we merge them
        Collections.sort(mPotentialBreakEnds,
                Comparator.comparing(TelomericBreakEnd::getType)
                .thenComparing(TelomericBreakEnd::getChromosome, ContigComparator.INSTANCE)
                .thenComparingInt(TelomericBreakEnd::getPosition));

        // now we can go through this list and merge the ones that are the same or similar
        CandidateBreakEnd currentBreakEnd = null;

        for (TelomericBreakEnd tbe : mPotentialBreakEnds)
        {
            if (currentBreakEnd != null &&
                currentBreakEnd.Type == tbe.getType() &&
                currentBreakEnd.Chromosome.equals(tbe.getChromosome()) &&
                Math.abs(currentBreakEnd.FirstPosition - tbe.getPosition()) <= mFuzzyMatchDistance)
            {
                // if they have same type and same chromosome, we check if they are nearby
                // we can combine them
                currentBreakEnd.LastPosition = tbe.getPosition();
                currentBreakEnd.NumSupportingFragments++;
            }
            else
            {
                // make a new one
                currentBreakEnd = new CandidateBreakEnd();
                currentBreakEnd.Type = tbe.getType();
                currentBreakEnd.Chromosome = tbe.getChromosome();
                currentBreakEnd.FirstPosition = tbe.getPosition();
                currentBreakEnd.FirstPosition = tbe.getPosition();
                consolidatedBreakEnds.add(currentBreakEnd);
            }
        }

        // after we do this, we want to check that we have sufficient evidence that this break end is real
        // what we want to do is to maybe go through all the reads again and find if any fragments that
        // span this break contains one side telomere and the other side non telomere

        for (CandidateBreakEnd candidateBreakEnd: consolidatedBreakEnds)
        {
            if (candidateBreakEnd.NumSupportingFragments >= mMinSupportingFragments)
            {
                int midPosition = candidateBreakEnd.FirstPosition + (candidateBreakEnd.LastPosition - candidateBreakEnd.FirstPosition) / 2;
                TelomericBreakEnd tbe = new TelomericBreakEnd(candidateBreakEnd.Type, candidateBreakEnd.Chromosome, midPosition);
                TE_LOGGER.trace("telomeric break end: {}, num supporting fragments: {}", tbe, candidateBreakEnd.NumSupportingFragments);
                mConsolidatedBreakEnds.add(tbe);
            }
        }

        GenomePosition[] positions = {
                GenomePositions.create("2", 241453197),
                GenomePositions.create("2", 242944598),
                GenomePositions.create("3", 88789845),
                GenomePositions.create("4", 127623126),
                GenomePositions.create("5", 149416338),
                GenomePositions.create("6", 16683178),
                GenomePositions.create("6", 169781199),
                GenomePositions.create("17", 22208891),
                GenomePositions.create("19", 53780521),
                GenomePositions.create("22", 21384627),
                GenomePositions.create("X", 12235971),
                GenomePositions.create("X", 19413088)
        };

        // check if any of them are missing
        for (GenomePosition p : positions)
        {
            boolean found = false;
            for (CandidateBreakEnd candidateBreakEnd: consolidatedBreakEnds)
            {
                if (candidateBreakEnd.FirstPosition == p.position() &&
                candidateBreakEnd.Chromosome.equals(p.chromosome()))
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                TE_LOGGER.info("missing: {}", p);
            }
        }
    }

    // look at soft clip site, and one part is mapped to a genome, and the
    // other part is telomeric

    // A right soft clip looks like 73S78M
    // A left soft clip looks like 78M32S

    // we want to work out whether one side is telomeric and the other side is not
    public TelomericBreakEnd findTelomericBreakEnd(@NotNull final SAMRecord r)
    {
        // TE_LOGGER.trace("record: {}", r);

        if (isInExcludedBaseRegion(r.getReferenceName(), r.getAlignmentStart(), r.getAlignmentEnd()))
        {
            return null;
        }

        String readString = r.getReadString();

        int leftClip = SamRecordUtils.leftSoftClip(r);
        int rightClip = SamRecordUtils.rightSoftClip(r);

        // now we work out if the clipped part is telomere
        String leftClipBases = readString.substring(0, leftClip);
        String rightClipBases = readString.substring(readString.length() - rightClip);

        String alignedBases = readString.substring(leftClip, readString.length() - rightClip);

        //String

        assert(leftClipBases.length() == leftClip);
        assert(rightClipBases.length() == rightClip);
        assert(leftClipBases.length() + rightClipBases.length() + alignedBases.length() == readString.length());

        if (!isStrictlyNonTelomeric(alignedBases) || TeloUtils.isPolyGC(alignedBases))
        {
            // aligned part must not be telomeric or poly G
            return null;
        }

        TelomericBreakEnd tbe = null;

        if (leftClipBases.length() >= mMinSoftClipLength + mBoundaryZone)
        {
            // we remove boundary zone
            String modLeftClipBases = leftClipBases.substring(0, leftClipBases.length() - mBoundaryZone);

            if (TelomereMatcher.calcGTelomereMatch(modLeftClipBases) >= mTelomereMatchThreshold)
            {
                tbe = new TelomericBreakEnd(TelomericBreakEnd.Type.LEFT_G_TELOMERIC, r.getReferenceName(),  r.getAlignmentStart());
            }
            else if (TelomereMatcher.calcCTelomereMatch(modLeftClipBases) >= mTelomereMatchThreshold)
            {
                tbe = new TelomericBreakEnd(TelomericBreakEnd.Type.LEFT_C_TELOMERIC, r.getReferenceName(), r.getAlignmentStart());
            }
        }
        if (rightClipBases.length() >= mMinSoftClipLength + mBoundaryZone)
        {
            // we remove boundary zone
            String modRightClipBases = rightClipBases.substring(mBoundaryZone);

            if (TelomereMatcher.calcGTelomereMatch(modRightClipBases) >= mTelomereMatchThreshold)
            {
                tbe = new TelomericBreakEnd(TelomericBreakEnd.Type.RIGHT_G_TELOMERIC, r.getReferenceName(), r.getAlignmentEnd());
            }
            else if (TelomereMatcher.calcCTelomereMatch(modRightClipBases) >= mTelomereMatchThreshold)
            {
                tbe = new TelomericBreakEnd(TelomericBreakEnd.Type.RIGHT_C_TELOMERIC, r.getReferenceName(), r.getAlignmentEnd());
            }
        }

        if (tbe != null)
        {
            TE_LOGGER.trace("record({}) -v strand({}) breakend({}) cigar({}) leftClip({}) rightClip({}) aligned({})",
                    r, r.getReadNegativeStrandFlag(), tbe, r.getCigarString(), leftClipBases, rightClipBases, alignedBases);
        }
        return tbe;
    }

    static boolean isStrictlyNonTelomeric(@NotNull final String seq)
    {
        return !seq.contains(TeloConstants.CANONICAL_TELOMERE_SEQ) && !seq.contains(TeloConstants.CANONICAL_TELOMERE_SEQ_REV);
    }

    // where 1 end maps in the POLY-G region of LINC00486 (v38: chr2:32,916,190-32,916,630; v37: 2:33,141,260-33,141,700).
    static boolean isInExcludedBaseRegion(String chromosome, int startPos, int endPos)
    {
        for (ChrBaseRegion excludedRegion : TeloConstants.EXCLUDED_BASE_REGIONS)
        {
            if (ContigComparator.INSTANCE.compare(chromosome, excludedRegion.chromosome()) == 0 &&
                    (excludedRegion.containsPosition(startPos) || excludedRegion.containsPosition(endPos)))
            {
                return true;
            }
        }

        return false;
    }
}
