package com.hartwig.hmftools.common.variant;

import org.apache.commons.math3.util.Pair;
import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

public final class VariantUtils
{
    public static Pair<Integer, String> relativePositionAndRef(final IndexedFastaSequenceFile reference, final VariantContext variant)
    {
        final int refLength = variant.getReference().getBaseString().length();
        final SAMSequenceRecord samSequenceRecord = reference.getSequenceDictionary().getSequence(variant.getContig());
        if(samSequenceRecord == null)
        {
            return new Pair<>(-1, Strings.EMPTY);
        }

        final int chromosomeLength = samSequenceRecord.getSequenceLength();
        int positionBeforeEvent = variant.getStart();

        int start = Math.max(positionBeforeEvent - 100, 1);
        int end = Math.min(positionBeforeEvent + refLength + 100 - 1, chromosomeLength - 1);
        int relativePosition = positionBeforeEvent - start;
        final String sequence;
        if(start < chromosomeLength && end < chromosomeLength)
        {
            sequence = reference.getSubsequenceAt(variant.getContig(), start, end).getBaseString();
        }
        else
        {
            sequence = Strings.EMPTY;
        }

        return new Pair<>(relativePosition, sequence);
    }

}
