package com.hartwig.hmftools.linx.neoepitope;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class RefGenomeSource implements RefGenomeInterface
{
    private final IndexedFastaSequenceFile mRefGenome;

    public RefGenomeSource(final IndexedFastaSequenceFile refGenome)
    {
        mRefGenome = refGenome;
    }

    @Override
    public String getBaseString(final String chromosome, long posStart, long posEnd)
    {
        return mRefGenome.getSubsequenceAt(chromosome, posStart, posEnd).getBaseString();

        /*
        final int refLength = variant.getReference().getBaseString().length();
        final int chromosomeLength = reference.getSequenceDictionary().getSequence(variant.getContig()).getSequenceLength();
        long positionBeforeEvent = variant.getStart();

        long start = Math.max(positionBeforeEvent - 100, 1);
        long end = Math.min(positionBeforeEvent + refLength + 100 - 1, chromosomeLength - 1);
        int relativePosition = (int) (positionBeforeEvent - start);
        if (start < chromosomeLength && end < chromosomeLength) {
            sequence =
        } else {
            sequence = Strings.EMPTY;
            LOGGER.warn("Requested base sequence outside of chromosome region!");
        }
        return new Pair<>(relativePosition, sequence);
        */
    }

}
