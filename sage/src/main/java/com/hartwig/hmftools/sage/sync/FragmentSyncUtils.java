package com.hartwig.hmftools.sage.sync;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;

public class FragmentSyncUtils
{
    protected static boolean isDeleteOrSplit(final CigarOperator element)
    {
        return element == D || element == N;
    }

    protected static boolean ignoreCigarOperatorMismatch(final CigarOperator first, final CigarOperator second)
    {
        return (first == M || first == S) && (second == M || second == S);
    }

    protected static boolean switchSoftClipToAligned(final CigarElement first, final CigarElement second)
    {
        if(first == null || second == null)
            return false;

        return (first.getOperator() == M && second.getOperator() == S) || (first.getOperator() == S && second.getOperator() == M);
    }

    protected static boolean overlappingCigarDiffs(final Cigar firstCigar, int firstPosStart, final Cigar secondCigar, int secondPosStart)
    {
        int firstAdjustedElementPosEnd = 0;
        int readPos = firstPosStart;
        for(CigarElement element : firstCigar.getCigarElements())
        {
            switch(element.getOperator())
            {
                case M:
                    readPos += element.getLength();
                    break;
                case D:
                case N:
                    readPos += element.getLength();
                    firstAdjustedElementPosEnd = readPos + element.getLength();
                    break;
                case I:
                    firstAdjustedElementPosEnd = readPos + 1;
                default:
                    break;
            }
        }

        int secondAdjustedElementPosStart = secondPosStart;
        for(CigarElement element : secondCigar.getCigarElements())
        {
            if(element.getOperator() == M)
                secondAdjustedElementPosStart += element.getLength();
            else if(element.getOperator().isIndel())
                break;
        }

        return firstAdjustedElementPosEnd >= secondAdjustedElementPosStart;
    }

    protected static boolean compatibleCigars(final Cigar firstCigar, final Cigar secondCigar)
    {
        // SC at the start and end are optional, but otherwise all elements must match length and type
        int j = 0;
        int i = 0;

        CigarElement firstElement = firstCigar.getCigarElements().get(i);
        CigarElement secondElement = secondCigar.getCigarElements().get(j);

        if(firstElement.getOperator() == S)
            ++i;

        if(secondElement.getOperator() == S)
            ++j;

        while(true)
        {
            firstElement = i < firstCigar.getCigarElements().size() ? firstCigar.getCigarElements().get(i) : null;
            secondElement = j < secondCigar.getCigarElements().size() ? secondCigar.getCigarElements().get(j) : null;

            if(firstElement == null && secondElement == null)
                break;

            if(firstElement == null)
                return secondElement.getOperator() == S;
            else if(secondElement == null)
                return firstElement.getOperator() == S;

            // must match types and lengths if not an alignment
            if(firstElement.getOperator() != secondElement.getOperator())
                return false;

            if(firstElement.getOperator() == S)
                return true;

            if(firstElement.getOperator() != M && firstElement.getLength() != secondElement.getLength())
                return false;

            ++i;
            ++j;
        }

        return true;
    }

    protected static byte[] getCombinedBaseAndQual(byte firstBase, byte firstQual, byte secondBase, byte secondQual)
    {
        if(firstBase == secondBase)
        {
            byte qual = (byte)max(firstQual, secondQual);
            return new byte[] { firstBase, qual };
        }
        else if(firstQual > secondQual)
        {
            // use the difference in quals
            return new byte[] { firstBase, (byte)((int)firstQual - (int)secondQual) };
        }
        else
        {
            return new byte[] { secondBase, (byte)((int)secondQual - (int)firstQual) };
        }
    }

    protected static SAMRecord buildSyncedRead(
            final SAMRecord referenceRead, final int combinedPosStart, final int combinedPosEnd,
            final byte[] combinedBases, final byte[] combinedBaseQualities,
            final List<CigarElement> combinedCigar, final int[] trucatedBases)
    {
        if(trucatedBases == null || (trucatedBases[0] == 0 && trucatedBases[0] == 0))
        {
            return buildSyncedRead(referenceRead, combinedPosStart, combinedBases, combinedBaseQualities, combinedCigar);
        }

        // bring in the read bases, quals, coords and cigar to the truncated positions
        int totalOffset = trucatedBases[0] + trucatedBases[1];
        int truncLength = combinedBases.length - totalOffset;
        byte[] truncBases = new byte[truncLength];
        byte[] truncBaseQuals = new byte[truncLength];

        for(int i = 0; i < truncLength; ++i)
        {
            truncBases[i] = combinedBases[i + trucatedBases[0]];
            truncBaseQuals[i] = combinedBaseQualities[i + trucatedBases[0]];
        }

        int truncatedFragmentStart = trucatedBases[2];

        for(int i = 0; i <= 1; ++i)
        {
            int offset = trucatedBases[i];

            while(offset > 0)
            {
                int cigarIndex = (i == 0) ? 0 : combinedCigar.size() - 1;
                CigarElement currentElement = combinedCigar.get(cigarIndex);

                if(offset < currentElement.getLength())
                {
                    combinedCigar.set(cigarIndex, new CigarElement(
                            currentElement.getLength() - offset, currentElement.getOperator()));
                    offset = 0;
                }
                else
                {
                    offset -= currentElement.getLength();
                    combinedCigar.remove(cigarIndex);
                }
            }
        }

        return buildSyncedRead(referenceRead, truncatedFragmentStart, truncBases, truncBaseQuals, combinedCigar);
    }

    protected static SAMRecord buildSyncedRead(
            final SAMRecord referenceRead, final int combinedPosStart, final byte[] combinedBases, final byte[] combinedBaseQualities,
            final List<CigarElement> combinedCigar)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        SAMRecord combinedRecord = recordBuilder.addFrag(
                referenceRead.getReadName(),
                referenceRead.getReferenceIndex(),
                combinedPosStart,
                referenceRead.getReadNegativeStrandFlag(),
                false,
                new Cigar(combinedCigar).toString(), "", 1, false);

        combinedRecord.setReadBases(combinedBases);
        combinedRecord.setAlignmentStart(combinedPosStart);
        combinedRecord.setReferenceIndex(referenceRead.getReferenceIndex());

        combinedRecord.setBaseQualities(combinedBaseQualities);
        combinedRecord.setReferenceName(referenceRead.getReferenceName());
        combinedRecord.setMateAlignmentStart(combinedPosStart);
        combinedRecord.setMateReferenceName(referenceRead.getReferenceName());
        combinedRecord.setMateReferenceIndex(referenceRead.getReferenceIndex());

        combinedRecord.setFlags(referenceRead.getFlags());

        combinedRecord.setMappingQuality(referenceRead.getMappingQuality());

        // no need to compute since both records have the same value and it remains unchanged
        combinedRecord.setInferredInsertSize(abs(referenceRead.getInferredInsertSize()));

        for(SAMRecord.SAMTagAndValue tagAndValue : referenceRead.getAttributes())
        {
            combinedRecord.setAttribute(tagAndValue.tag, tagAndValue.value);
        }

        return combinedRecord;
    }
}
