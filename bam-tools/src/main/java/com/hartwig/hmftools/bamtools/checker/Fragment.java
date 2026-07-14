package com.hartwig.hmftools.bamtools.checker;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_INDEX;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CHROMOSOME_NAME;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NO_CIGAR;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.READ_GROUP_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.ALIGNMENTS_DELIM;

import static htsjdk.samtools.CigarOperator.H;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamReadLite;
import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.Arrays;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

public class Fragment
{
    private final String mReadId;
    private List<SAMRecord> mReads;

    private String mReadGroupIdId;
    private List<BamReadLite> mLiteReads;

    // counts of supplementary and secondary reads
    private short mExpectedSupplementaryCount;
    private short mReceivedSupplementaryCount;

    private String mFirstPrimaryCigar;
    private String mSecondPrimaryCigar;

    // cached read info on primaries for supplementary hard-clip conversion
    private byte[] mFirstReadBases;
    private byte[] mSecondReadBases;
    private byte[] mFirstReadQuals;
    private byte[] mSecondReadQuals;
    private boolean mFirstNegOrientation;
    private boolean mSecondNegOrientation;

    private boolean mFirstUnmapped;
    private boolean mSecondUnmapped;

    private boolean mMateCigarFixed;

    public Fragment(final SAMRecord read)
    {
        mReadId = read.getReadName();
        mReads = Lists.newArrayListWithExpectedSize(2);
        mLiteReads = null;
        mExpectedSupplementaryCount = 0;
        mReceivedSupplementaryCount = 0;
        mFirstPrimaryCigar = null;
        mSecondPrimaryCigar = null;
        mMateCigarFixed = false;
        mFirstUnmapped = false;
        mSecondUnmapped = false;

        mFirstReadBases = null;
        mSecondReadBases = null;
        mFirstReadQuals = null;
        mSecondReadQuals = null;

        addRead(read);
    }

    public String readId() { return mReadId; }
    public List<SAMRecord> reads() { return mReads != null ? mReads : Collections.emptyList(); }
    public int readCount() { return mReads != null ? mReads.size() : (mLiteReads != null ? mLiteReads.size() : 0); }
    public short expectedSupplementaryCount() { return mExpectedSupplementaryCount; }
    public short receivedSupplementaryCount() { return mReceivedSupplementaryCount; }

    public String firstPrimaryCigar() { return mFirstPrimaryCigar; }
    public String secondPrimaryCigar() { return mSecondPrimaryCigar; }

    public byte[] firstReadBases() { return mFirstReadBases; }
    public byte[] firstBaseQuals() { return mFirstReadQuals; }
    public boolean firstNegOrientation() { return mFirstNegOrientation; }
    public byte[] secondReadBases() { return mSecondReadBases; }
    public byte[] secondBaseQuals() { return mSecondReadQuals; }
    public boolean secondNegOrientation() { return mSecondNegOrientation; }

    public boolean firstUnmapped() { return mFirstUnmapped; }
    public boolean secondUnmapped() { return mSecondUnmapped; }
    public int unmappedPrimaryCount() { return (mFirstUnmapped ? 1 : 0) + (mSecondUnmapped ? 1 : 0); }

    public void addRead(final SAMRecord read)
    {
        if(mReads == null)
            mReads = Lists.newArrayListWithExpectedSize(2);

        if(read.getSupplementaryAlignmentFlag())
        {
            if((mFirstUnmapped && read.getFirstOfPairFlag()) || (mSecondUnmapped && !read.getFirstOfPairFlag()))
                return;

            ++mReceivedSupplementaryCount;
        }
        else
        {
            if(read.getFirstOfPairFlag())
            {
                mFirstUnmapped = belowMinAlignmentScore(read);
                mFirstPrimaryCigar = !mFirstUnmapped ? read.getCigarString() : NO_CIGAR;

            }
            else
            {
                mSecondUnmapped = belowMinAlignmentScore(read);
                mSecondPrimaryCigar = !mSecondUnmapped ? read.getCigarString() : NO_CIGAR;
            }

            if(mFirstUnmapped || mSecondUnmapped)
                purgeSupplementaryReads();

            checkSupplementaryData(read);
        }

        mReads.add(read);
    }

    private boolean belowMinAlignmentScore(final SAMRecord read)
    {
        Integer asScore = read.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
        return asScore != null && asScore.intValue() < CheckConfig.Params.MinAlignmentScore;
    }

    private void checkSupplementaryData(final SAMRecord read)
    {
        if((mFirstUnmapped && read.getFirstOfPairFlag()) || (mSecondUnmapped && !read.getFirstOfPairFlag()))
        {
            // clear supp attributes
            read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
            return;
        }

        // filter out supps with low alignment scores, rebuilding the SA tag if required
        List<SupplementaryReadData> suppDataList = SupplementaryReadData.extractAlignments(read);

        if(suppDataList == null)
            return;

        List<SupplementaryReadData> validSuppDataList = Lists.newArrayListWithExpectedSize(suppDataList.size());

        for(SupplementaryReadData suppData : suppDataList)
        {
            Cigar suppCigar = CigarUtils.cigarFromStr(suppData.Cigar);
            int alignedBases = suppCigar.getCigarElements().stream()
                    .filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();

            if(alignedBases >= CheckConfig.Params.MinAlignmentScore)
            {
                validSuppDataList.add(suppData);
            }
        }

        if(validSuppDataList.isEmpty())
        {
            read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, null);
            return;
        }

        mExpectedSupplementaryCount += validSuppDataList.size();

        if(validSuppDataList.size() < suppDataList.size())
        {
            StringJoiner sj = new StringJoiner(ALIGNMENTS_DELIM);

            for(SupplementaryReadData suppData : validSuppDataList)
            {
                sj.add(suppData.asSamTag());
            }

            read.setAttribute(SUPPLEMENTARY_ATTRIBUTE, sj.toString());
        }
    }

    public boolean hasPrimaryInfo() { return mFirstPrimaryCigar != null && mSecondPrimaryCigar != null; }
    public boolean requiredMateCigarFix() { return mMateCigarFixed; }

    public boolean isComplete()
    {
        return hasPrimaryInfo() && mExpectedSupplementaryCount == mReceivedSupplementaryCount;
    }

    public synchronized boolean mergeFragment(
            final Fragment fragment, final SAMFileHeader samFileHeader, final List<SAMRecord> completeReads)
    {
        // copy reads or primary cigar info if the fragment has already written the primaries
        transfer(fragment);

        // keep the fragment's reads if one or both primaries are missing
        if(!hasPrimaryInfo())
        {
            serialiseReads();
            return false;
        }

        // primary cigar info is complete, so write any cached reads
        deserialiseReads(samFileHeader);

        List<SAMRecord> fragCompleteReads = extractCompleteReads();

        completeReads.addAll(fragCompleteReads);

        // keep just its primary info if its is waiting on supplementaries only
        if(isComplete())
            return true;

        serialiseReads();
        return false;
    }

    public void transfer(final Fragment other)
    {
        other.reads().forEach(x -> addRead(x));

        if(other.hasPrimaryInfo() && !hasPrimaryInfo())
        {
            mFirstPrimaryCigar = other.firstPrimaryCigar();
            mSecondPrimaryCigar = other.secondPrimaryCigar();
            mFirstUnmapped = other.firstUnmapped();
            mSecondUnmapped = other.secondUnmapped();
            mExpectedSupplementaryCount += other.expectedSupplementaryCount();
        }

        if(other.firstReadBases() != null && mFirstReadBases == null)
        {
            mFirstReadBases = other.firstReadBases();
            mFirstReadQuals = other.firstBaseQuals();
            mFirstNegOrientation = other.firstNegOrientation();
        }

        if(other.secondReadBases() != null && mSecondReadBases == null)
        {
            mSecondReadBases = other.secondReadBases();
            mSecondReadQuals = other.secondBaseQuals();
            mSecondNegOrientation = other.secondNegOrientation();
        }
    }

    public List<SAMRecord> extractCompleteReads()
    {
        if(!hasPrimaryInfo() || mReads == null || mReads.isEmpty())
            return Collections.emptyList();

        checkPrimaryUnmapping();

        setMateCigar();

        convertHardClips();

        List<SAMRecord> reads = Lists.newArrayList(mReads);
        mReads.clear();
        return reads;
    }

    private void setMateCigar()
    {
        for(SAMRecord read : mReads)
        {
            mMateCigarFixed |= !read.hasAttribute(MATE_CIGAR_ATTRIBUTE);

            if(read.getFirstOfPairFlag())
            {
                read.setAttribute(MATE_CIGAR_ATTRIBUTE, mSecondPrimaryCigar);
            }
            else
            {
                read.setAttribute(MATE_CIGAR_ATTRIBUTE, mFirstPrimaryCigar);
            }
        }
    }

    public void cachePrimaryBaseInfo(final SAMRecord read)
    {
        if(read.getSupplementaryAlignmentFlag() || !read.hasAttribute(SUPPLEMENTARY_ATTRIBUTE))
            return;

        if(read.getFirstOfPairFlag())
        {
            mFirstReadBases = read.getReadBases();
            mFirstReadQuals = read.getBaseQualities();
            mFirstNegOrientation = read.getReadNegativeStrandFlag();
        }
        else
        {
            mSecondReadBases = read.getReadBases();
            mSecondReadQuals = read.getBaseQualities();
            mSecondNegOrientation = read.getReadNegativeStrandFlag();
        }
    }

    private void purgeSupplementaryReads()
    {
        int index = 0;

        while(index < mReads.size())
        {
            SAMRecord read = mReads.get(index);

            if(mFirstUnmapped && read.getFirstOfPairFlag() && read.getSupplementaryAlignmentFlag())
            {
                --mReceivedSupplementaryCount;
                mReads.remove(index);
                continue;
            }
            else if(mSecondUnmapped && !read.getFirstOfPairFlag() && read.getSupplementaryAlignmentFlag())
            {
                --mReceivedSupplementaryCount;
                mReads.remove(index);
                continue;
            }

            ++index;
        }
    }

    private void checkPrimaryUnmapping()
    {
        if(!mFirstUnmapped && !mSecondUnmapped)
            return;

        SAMRecord firstPrimary = null;
        SAMRecord secondPrimary = null;

        for(SAMRecord read : mReads)
        {
            if(read.getSupplementaryAlignmentFlag())
                continue;

            if(read.getFirstOfPairFlag())
            {
                firstPrimary = read;
            }
            else
            {
                secondPrimary = read;
            }
        }

        if(mFirstUnmapped && mSecondUnmapped)
        {
            for(SAMRecord read : mReads)
            {
                read.setProperPairFlag(false);
                read.setInferredInsertSize(0);
                read.setReadUnmappedFlag(true);
                read.setMateUnmappedFlag(true);
                read.setMappingQuality(0);
                read.setAlignmentStart(0);

                read.setReferenceIndex(NO_CHROMOSOME_INDEX);
                read.setReferenceName(NO_CHROMOSOME_NAME);
                read.setCigarString(NO_CIGAR);

                read.setMateAlignmentStart(0);
                read.setMateReferenceIndex(NO_CHROMOSOME_INDEX);
                read.setMateReferenceName(NO_CHROMOSOME_NAME);
                read.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
            }

            return;
        }

        SAMRecord mappedRead = !mFirstUnmapped ? firstPrimary : secondPrimary;
        SAMRecord unmappedRead = mFirstUnmapped ? firstPrimary : secondPrimary;

        if(mappedRead != null)
        {
            mappedRead.setProperPairFlag(false);

            mappedRead.setInferredInsertSize(0);
            mappedRead.setMateUnmappedFlag(true);
            mappedRead.setMateAlignmentStart(mappedRead.getAlignmentStart());
            mappedRead.setMateReferenceIndex(mappedRead.getReferenceIndex());
            mappedRead.setAttribute(MATE_CIGAR_ATTRIBUTE, null);
        }

        if(unmappedRead != null)
        {
            unmappedRead.setProperPairFlag(false);
            unmappedRead.setInferredInsertSize(0);
            unmappedRead.setReadUnmappedFlag(true);
            unmappedRead.setMappingQuality(0);
            unmappedRead.setReferenceIndex(mappedRead.getReferenceIndex());
            unmappedRead.setAlignmentStart(mappedRead.getAlignmentStart());
            unmappedRead.setCigarString(NO_CIGAR);
        }
    }

    private void convertHardClips()
    {
        if(mFirstReadBases == null && mSecondReadBases == null)
            return;

        for(SAMRecord read : mReads)
        {
            if(!read.getSupplementaryAlignmentFlag())
                continue;

            if(read.getCigar().getCigarElements().stream().noneMatch(x -> x.getOperator() == H))
                continue;

            if(read.getFirstOfPairFlag())
            {
                convertHardClips(read, mFirstReadBases, mFirstReadQuals, mFirstNegOrientation);
            }
            else
            {
                convertHardClips(read, mSecondReadBases, mSecondReadQuals, mSecondNegOrientation);
            }
        }
    }

    @VisibleForTesting
    public static void convertHardClips(final SAMRecord read, final byte[] readBases, final byte[] baseQuals, final boolean negOrientation)
    {
        List<CigarElement> cigarElements = read.getCigar().getCigarElements();
        List<CigarElement> newCigarElements = Lists.newArrayListWithExpectedSize(cigarElements.size());

        if(negOrientation != read.getReadNegativeStrandFlag())
        {
            read.setReadBases(Nucleotides.reverseComplementBases(readBases));
            read.setBaseQualities(Arrays.reverseArray(baseQuals));
        }
        else
        {
            read.setReadBases(readBases);
            read.setBaseQualities(baseQuals);
        }

        for(CigarElement cigarElement : cigarElements)
        {
            CigarOperator operator = cigarElement.getOperator() == H ? S : cigarElement.getOperator();
            newCigarElements.add(new CigarElement(cigarElement.getLength(), operator));
        }

        read.setCigar(new Cigar(newCigarElements));
    }

    public String toString()
    {
        return format("id(%s) reads(%d) primary(first=%s second=%s) supps(expected=%d received=%d)",
                mReadId, readCount(), mFirstPrimaryCigar != null, mSecondPrimaryCigar != null,
                mExpectedSupplementaryCount, mReceivedSupplementaryCount);
    }

    public synchronized void serialiseReads()
    {
        if(mReads == null || mReads.isEmpty())
            return;

        if(mLiteReads == null)
            mLiteReads = Lists.newArrayListWithExpectedSize(mReads.size());

        if(mReadGroupIdId == null)
            mReadGroupIdId = mReads.get(0).getStringAttribute(READ_GROUP_ATTRIBUTE);

        for(SAMRecord read : mReads)
        {
            mLiteReads.add(new BamReadLite(read, true));
        }

        mReads.clear();
        mReads = null;
    }

    public void deserialiseReads(final SAMFileHeader samFileHeader)
    {
        if(mLiteReads == null || mLiteReads.isEmpty())
            return;

        if(mReads == null)
            mReads = Lists.newArrayListWithExpectedSize(mLiteReads.size());

        for(BamReadLite readLite : mLiteReads)
        {
            SAMRecord read = BamReadLite.from(readLite, samFileHeader, mReadId, mReadGroupIdId);
            mReads.add(read);
        }

        mLiteReads.clear();
        mLiteReads = null;
    }
}
