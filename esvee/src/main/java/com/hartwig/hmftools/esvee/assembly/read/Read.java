package com.hartwig.hmftools.esvee.assembly.read;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.cigarElementsToStr;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_POSITION;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getMateAlignmentEnd;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.extractAlignment;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.utils.Arrays.copyArray;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.READ_ID_TRIMMER;
import static com.hartwig.hmftools.esvee.assembly.read.ReadAdjustments.LOW_QUAL_SCORE;
import static com.hartwig.hmftools.esvee.common.CommonUtils.belowMinQual;
import static com.hartwig.hmftools.esvee.common.IndelCoords.findIndelCoords;
import static com.hartwig.hmftools.esvee.common.SvConstants.BAM_HEADER_SAMPLE_INDEX_TAG;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_SUPPORT_LENGTH;

import static htsjdk.samtools.CigarOperator.S;
import static htsjdk.samtools.util.StringUtil.bytesToString;

import java.util.List;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.common.IndelCoords;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class Read
{
    private final SAMRecord mRecord;

    // cached state and adjusted properties of the read
    private final String mOrigCigarString;
    private String mCigarString;
    private List<CigarElement> mCigarElements;

    private String mId; // initialised and trimmed on demand
    private int mAlignmentStart;
    private int mAlignmentEnd;
    private int mUnclippedStart;
    private int mUnclippedEnd;
    private Integer mMateAlignmentEnd;
    private byte[] mBases;
    private byte[] mBaseQuals;

    // fragment state
    private Read mMateRead;
    private boolean mHasJunctionMate; // mate read is a split/junction read
    private boolean mSuppDataExtracted;
    private SupplementaryReadData mSupplementaryData;

    private boolean mCheckedIndelCoords;

    private IndelCoords mIndelCoords;
    private boolean mInvalidIndel;

    // for reads with indels >= min-length, the inferred unclipped adjustment from converting them all
    private Integer mIndelInferredUnclippedStart;
    private Integer mIndelInferredUnclippedEnd;

    private boolean mIsReference;
    private boolean mHasLineTail;
    private int mTrimCount; // if non-zero this has been done to the 3' end
    private boolean mLowQualTrimmed;

    public Read(final SAMRecord record)
    {
        mRecord = record;

        mId = null;

        mOrigCigarString = record.getCigarString();
        mCigarString = null;
        mCigarElements = cigarElementsFromStr(mOrigCigarString);

        setBoundaries(mRecord.getAlignmentStart());
        mBases = null;
        mBaseQuals = null;
        mMateAlignmentEnd = null;
        mIsReference = false;
        mMateRead = null;
        mHasJunctionMate = false;
        mSuppDataExtracted = false;
        mSupplementaryData = null;
        mCheckedIndelCoords = false;

        // only set for adjusted indel reads
        mIndelCoords = null;
        mIndelInferredUnclippedStart = null;
        mIndelInferredUnclippedEnd = null;
        mInvalidIndel = false;

        mHasLineTail = false;
        mTrimCount = 0;
        mLowQualTrimmed = false;
    }

    private void setBoundaries(int newReadStart)
    {
        mAlignmentStart = newReadStart;
        mUnclippedStart = mAlignmentStart;

        if(mCigarElements.isEmpty())
        {
            // undefined for unmapped reads
            mAlignmentEnd = mAlignmentStart;
            mUnclippedEnd = mAlignmentStart;
            return;
        }

        int currentPosition = mAlignmentStart;

        for(int i = 0; i < mCigarElements.size(); ++i)
        {
            CigarElement element = mCigarElements.get(i);

            if(i == 0 && element.getOperator() == S)
                mUnclippedStart -= element.getLength();

            if(element.getOperator().consumesReferenceBases())
                currentPosition += element.getLength();

            if(i == mCigarElements.size() - 1)
            {
                mAlignmentEnd = currentPosition - 1;
                mUnclippedEnd = element.getOperator() == S ? mAlignmentEnd + element.getLength() : mAlignmentEnd;
            }
        }
    }

    public SAMRecord bamRecord() { return mRecord; }

    public void setMateRead(final Read mate)
    {
        mMateRead = mate;
        mMateAlignmentEnd = mate.alignmentEnd();
    }

    public Read mateRead() { return mMateRead; }

    public boolean hasJunctionMate() { return mHasJunctionMate; }
    public void markJunctionMate() { mHasJunctionMate = true; }

    public String id()
    {
        if(mId == null)
            mId = READ_ID_TRIMMER.trim(mRecord.getReadName());

        return mId;
    }

    public String chromosome() { return mRecord.getReferenceName(); }

    public List<CigarElement> cigarElements() { return mCigarElements; }
    public String cigarString() { return mCigarString != null ? mCigarString : mOrigCigarString; }
    public String originalCigarString() { return mOrigCigarString; }
    private void updateCigarString() { mCigarString = cigarElementsToStr(mCigarElements); }

    public int alignmentStart() { return mAlignmentStart; }
    public int alignmentEnd() { return mAlignmentEnd; }

    public int unclippedStart()  { return mUnclippedStart; }
    public int unclippedEnd() { return mUnclippedEnd; }

    public int fragmentEnd() { return positiveStrand() ? mUnclippedStart : mUnclippedEnd; }

    // convenience
    public boolean isLeftClipped() { return mUnclippedStart != mAlignmentStart || mIndelInferredUnclippedStart != null; }
    public boolean isRightClipped() { return mUnclippedEnd != mAlignmentEnd || mIndelInferredUnclippedEnd != null; }

    public int leftClipLength() { return max(mAlignmentStart - mUnclippedStart, 0); } // no known need to use the indel-implied SC value
    public int rightClipLength() { return max(mUnclippedEnd - mAlignmentEnd, 0); }

    public byte[] getBases() { return mBases != null ? mBases : mRecord.getReadBases(); }
    public byte[] getBaseQuality() { return mBaseQuals != null ? mBaseQuals : mRecord.getBaseQualities(); }
    public int basesLength() { return mBases != null ? mBases.length : mRecord.getReadBases().length; }

    // flags
    public int getFlags() { return mRecord.getFlags(); }
    public boolean isUnmapped() { return mRecord.getReadUnmappedFlag(); }
    public boolean isPairedRead() { return mRecord.getReadPairedFlag(); }

    public boolean positiveStrand() { return !mRecord.getReadNegativeStrandFlag(); }
    public boolean negativeStrand() { return mRecord.getReadNegativeStrandFlag(); }
    public Orientation orientation() { return mRecord.getReadNegativeStrandFlag() ? REVERSE : FORWARD; }

    public boolean firstInPair() { return mRecord.getReadPairedFlag() && mRecord.getFirstOfPairFlag(); }

    public int mappingQuality() { return mRecord.getMappingQuality(); }

    // if the mate read reference is set, use this so it can be clear on any unmapping state
    public String mateChromosome()
    {
        if(mMateRead != null)
            return mMateRead.chromosome();;

        return isMateMapped() ? mRecord.getMateReferenceName() : null;
    }

    public int mateAlignmentStart()
    {
        if(mMateRead != null)
            return mMateRead.alignmentStart();

        return mRecord.getMateAlignmentStart();
    }

    public int mateAlignmentEnd()
    {
        if(mMateRead != null)
            return mMateRead.alignmentEnd();

        if(mMateAlignmentEnd != null)
            return mMateAlignmentEnd;

        mMateAlignmentEnd = getMateAlignmentEnd(mRecord);

        if(mMateAlignmentEnd == NO_POSITION)
            mMateAlignmentEnd = mateAlignmentStart() + basesLength() - 1;

        return mMateAlignmentEnd;
    }

    public int mateFragmentEnd()
    {
        if(mMateRead != null)
            return mMateRead.orientation().isForward() ? mMateRead.unclippedStart() : mMateRead.unclippedEnd();

        String mateCigarStr = mRecord.getStringAttribute(MATE_CIGAR_ATTRIBUTE);

        if(mateCigarStr == null)
            return NO_POSITION;

        return getFivePrimeUnclippedPosition(mateAlignmentStart(), mateCigarStr, !mRecord.getMateNegativeStrandFlag());
    }

    public boolean isMateMapped()
    {
        if(mMateRead != null)
            return !mMateRead.isUnmapped();

        return mRecord.getReadPairedFlag() && !mRecord.getMateUnmappedFlag();
    }

    public boolean isMateUnmapped() { return !isMateMapped(); }

    public Orientation mateOrientation()
    {
        if(!mRecord.getReadPairedFlag())
            return FORWARD;

        return mRecord.getMateNegativeStrandFlag() ? REVERSE : FORWARD; }

    public boolean hasSupplementary() { return supplementaryData() != null; }
    public boolean isSupplementary() { return mRecord.getSupplementaryAlignmentFlag(); }

    public void makeReadLinks(final Read other)
    {
        if(mRecord.getSupplementaryAlignmentFlag() == other.bamRecord().getSupplementaryAlignmentFlag()
        && firstInPair() != other.firstInPair())
        {
            mMateRead = other;
            other.setMateRead(this);

            if(!bamRecord().getReadPairedFlag()) // deficiency in Redux, solved in v1.0
                bamRecord().setReadPairedFlag(true);

            if(!other.bamRecord().getReadPairedFlag())
                other.bamRecord().setReadPairedFlag(true);

            // correct unmapped flags if applicable - due to the lack of mate cigar impacting Redux
            if(isUnmapped() && !other.bamRecord().getMateUnmappedFlag())
                other.bamRecord().setMateUnmappedFlag(true);
            else if(other.isUnmapped() && !bamRecord().getMateUnmappedFlag())
                bamRecord().setMateUnmappedFlag(true);
        }
    }

    public SupplementaryReadData supplementaryData()
    {
        if(!mSuppDataExtracted)
        {
            mSuppDataExtracted = true;
            mSupplementaryData = extractAlignment(mRecord);
        }

        return mSupplementaryData;
    }

    public int getReadIndexAtReferencePosition(final int refPosition)
    {
        return getReadIndexAtReferencePosition(refPosition, false);
    }

    public int getReadIndexAtReferencePosition(final int refPosition, boolean allowExtrapolation)
    {
        return ReadUtils.getReadIndexAtReferencePosition(this, refPosition, allowExtrapolation);
    }

    public IndelCoords indelCoords()
    {
        if(!mCheckedIndelCoords && !mInvalidIndel)
        {
            mCheckedIndelCoords = true;
            mIndelCoords = findIndelCoords(mAlignmentStart, mCigarElements, MIN_INDEL_SUPPORT_LENGTH);
        }

        return mIndelCoords;
    }

    public void markInvalidIndel() { mInvalidIndel = true; }
    public boolean invalidIndel() { return mInvalidIndel; }

    public boolean matchesFragment(final Read other, boolean allowReadMatch)
    {
        if(!id().equals(other.id()))
            return false;

        return allowReadMatch || getFlags() != other.getFlags();
    }

    public static List<Read> findMatchingFragmentSupport(final List<Read> support, final Read read)
    {
        return support.stream().filter(x -> x.matchesFragment(read, false)).collect(Collectors.toList());
    }

    public String toString()
    {
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                id(), chromosome(), mAlignmentStart, mAlignmentEnd, cigarString(),
                mateChromosome(), mateAlignmentStart(), mRecord.getFlags());
    }

    public int sampleIndex()
    {
        String sampleIndex = mRecord.getHeader().getAttribute(BAM_HEADER_SAMPLE_INDEX_TAG);
        return sampleIndex != null ? Integer.parseInt(sampleIndex) : 0;
    }

    public boolean isReference() { return mIsReference; }
    public void markReference() { mIsReference = true; }

    public boolean hasLineTail() { return mHasLineTail; }
    public void markLineTail() { mHasLineTail = true; }


    @VisibleForTesting
    public String getBasesString() { return bytesToString(getBases()); }

    public void trimBases(int count, boolean fromStart)
    {
        int remainingBases = count;
        int newBaseLength = max(basesLength() - count, 1);
        byte[] newBases = new byte[newBaseLength];
        byte[] newBaseQuals = new byte[newBaseLength];

        int newReadStart = mAlignmentStart;

        if(fromStart)
        {
            while(mCigarElements.size() > 0)
            {
                CigarElement element = mCigarElements.get(0);

                if(remainingBases == 0)
                {
                    // remove non-read elements
                    if(element.getOperator().consumesReferenceBases() && !element.getOperator().consumesReadBases())
                    {
                        mCigarElements.remove(0);
                        newReadStart += element.getLength();
                        continue;
                    }

                    break;
                }

                if(!element.getOperator().consumesReadBases())
                {
                    mCigarElements.remove(0);

                    if(element.getOperator().consumesReferenceBases())
                        newReadStart += element.getLength();

                    continue;
                }

                if(element.getLength() <= remainingBases)
                {
                    mCigarElements.remove(0);
                    remainingBases -= element.getLength();

                    if(element.getOperator().consumesReferenceBases())
                        newReadStart += element.getLength();
                }
                else
                {
                    mCigarElements.set(0, new CigarElement(element.getLength() - remainingBases, element.getOperator()));

                    if(element.getOperator().consumesReferenceBases())
                        newReadStart += remainingBases;

                    remainingBases = 0;
                }
            }

            copyArray(getBases(), newBases, count, 0);
            copyArray(getBaseQuality(), newBaseQuals, count, 0);
        }
        else
        {
            while(mCigarElements.size() > 0)
            {
                int lastIndex = mCigarElements.size() - 1;
                CigarElement element = mCigarElements.get(lastIndex);

                if(remainingBases == 0)
                {
                    // remove non-read elements
                    if(element.getOperator().consumesReferenceBases() && !element.getOperator().consumesReadBases())
                    {
                        mCigarElements.remove(lastIndex);
                        continue;
                    }

                    break;
                }

                if(!element.getOperator().consumesReadBases())
                {
                    mCigarElements.remove(lastIndex);
                    continue;
                }

                if(element.getLength() <= remainingBases)
                {
                    mCigarElements.remove(lastIndex);
                    remainingBases -= element.getLength();
                }
                else
                {
                    mCigarElements.set(lastIndex, new CigarElement(element.getLength() - remainingBases, element.getOperator()));
                    remainingBases = 0;
                }
            }

            copyArray(getBases(), newBases, 0, 0);
            copyArray(getBaseQuality(), newBaseQuals, 0, 0);
        }

        mBases = newBases;
        mBaseQuals = newBaseQuals;
        mTrimCount = count;

        updateCigarString();
        setBoundaries(newReadStart);
    }

    public int baseTrimCount() { return mTrimCount; }

    public void markLowQualTrimmed() { mLowQualTrimmed = true; }
    public boolean lowQualTrimmed() { return mLowQualTrimmed; }

    public void setIndelInferredUnclippedPositions(@Nullable final Integer inferredStart, @Nullable final Integer inferredEnd)
    {
        mIndelInferredUnclippedStart = inferredStart;
        mIndelInferredUnclippedEnd = inferredEnd;
    }

    public boolean hasIndelImpliedUnclippedStart() { return mIndelInferredUnclippedStart != null; }
    public boolean hasIndelImpliedUnclippedEnd() { return mIndelInferredUnclippedEnd != null; }

    public int indelImpliedUnclippedStart() { return mIndelInferredUnclippedStart != null ? mIndelInferredUnclippedStart : 0; }
    public int indelImpliedUnclippedEnd() { return mIndelInferredUnclippedEnd != null ? mIndelInferredUnclippedEnd : 0; }

    // take indel implied read ends into consideration for methods requiring the maximum possible read soft-clip extension
    // note: converted INDELs from deletes may have their unclipped position inside the alignment
    public int minUnclippedStart()
    {
        return mIndelInferredUnclippedStart == null ? mUnclippedStart : min(mUnclippedStart, mIndelInferredUnclippedStart);
    }

    public int maxUnclippedEnd()
    {
        return mIndelInferredUnclippedEnd == null ? mUnclippedEnd : max(mUnclippedEnd, mIndelInferredUnclippedEnd);
    }

    public synchronized void trimLowQualBases()
    {
        if(lowQualTrimmed())
            return;

        boolean fromStart = negativeStrand();

        int baseLength = basesLength();
        int baseIndex = fromStart ? 0 : baseLength - 1;

        int checkedBases = 0;

        double lowestScore = 0;
        double currentScore = 0;
        int lastLowestScoreIndex = 0;

        while(baseIndex >= 0 && baseIndex < baseLength)
        {
            ++checkedBases;

            if(belowMinQual(getBaseQuality()[baseIndex]))
            {
                currentScore -= LOW_QUAL_SCORE;

                if(currentScore < lowestScore)
                    lastLowestScoreIndex = checkedBases;
            }
            else
            {
                ++currentScore;
            }

            if(fromStart)
                ++baseIndex;
            else
                --baseIndex;
        }

        if(lastLowestScoreIndex > 0)
        {
            trimBases(lastLowestScoreIndex, fromStart);
            markLowQualTrimmed();
        }
    }
}
