package com.hartwig.hmftools.sage.common;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.sage.SageConstants.SC_READ_EVENTS_FACTOR;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

public final class NumberEvents
{
    public static int calcAdjustedNumMutations(
            final SAMRecord record, final RefSequence refSequence, @Nullable final RefGenomeInterface refGenome)
    {
        Object nm = record.getAttribute(NUM_MUTATONS_ATTRIBUTE);

        if(nm != null && nm instanceof Integer)
            return (int)nm;

        return calcMutations(record, refSequence, refGenome);
    }

    public static int calcMutations(final SAMRecord record, final RefSequence refSequence, @Nullable final RefGenomeInterface refGenome)
    {
        String contig = record.getReferenceName();
        int errorCount = 0;
        byte[] readBases = record.getReadBases();

        for(AlignmentBlock alignment : record.getAlignmentBlocks())
        {
            int alignmentLength = alignment.getLength();
            int referenceStart = alignment.getReferenceStart();
            int referenceEnd = referenceStart + alignmentLength - 1;
            int readIndexStart = alignment.getReadStart() - 1; // since is 1-based

            byte[] refBases;

            if(referenceStart >= refSequence.Start && referenceEnd <= refSequence.End)
            {
                refBases = refSequence.baseRange(referenceStart, referenceEnd);
            }
            else
            {
                // can occur for RNA reads with N-splits spanning past the end of the ref sequence
                refBases = refGenome.getBases(contig, referenceStart, referenceEnd);
            }

            if(refBases == null || refBases.length != alignmentLength)
                break;

            int readIndex = readIndexStart;

            for(int refIndex = 0; refIndex < alignmentLength; ++refIndex, ++readIndex)
            {
                if(readIndex >= readBases.length)
                    break;

                if(readBases[readIndex] != refBases[refIndex])
                    ++errorCount;
            }
        }

        return errorCount;
    }

    public static double calcSoftClipAdjustment(final SAMRecord record)
    {
        return calcSoftClipAdjustment(leftSoftClipLength(record) + rightSoftClipLength(record));
    }

    public static double calcSoftClipAdjustment(int softClipLength)
    {
        return softClipLength > 0 ? max(1, softClipLength / SC_READ_EVENTS_FACTOR) : 0;
    }

    public static int calcWithMnvRaw(int numberOfEvents, final String ref, final String alt)
    {
        // number of events includes each SNV as an additional event - this unfairly penalises MNVs
        int differentBases = 0;
        for(int i = 0; i < alt.length(); i++)
        {
            if(alt.charAt(i) != ref.charAt(i))
            {
                differentBases++;
            }
        }

        // We subtract one later when we actually use this value so we need to add one back in here to be consistent with SNVs and INDELs
        return numberOfEvents - differentBases + 1;
    }
}
