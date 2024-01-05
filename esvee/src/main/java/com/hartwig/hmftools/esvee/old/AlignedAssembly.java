package com.hartwig.hmftools.esvee.old;

import static com.hartwig.hmftools.esvee.read.ReadUtils.flipRead;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.esvee.read.Read;

import htsjdk.samtools.util.SequenceUtil;

public class AlignedAssembly extends SupportedAssembly implements AlignedSequence
{
    public final GappedAssembly Source;
    private final List<Alignment> mAlignment;

    public AlignedAssembly(final GappedAssembly source, final List<Alignment> alignment)
    {
        this(source.Assembly, source, alignment);

        for(ReadSupport support : source.readSupport())
        {
            addEvidenceAt(support.Read, support.Index);
        }

        recalculateBaseQuality();
    }

    private AlignedAssembly(final String assembly, final GappedAssembly source, final List<Alignment> alignment)
    {
        super(source.Name, assembly);

        Source = source;
        mAlignment = alignment;
    }

    @Override
    public List<Alignment> getAlignmentBlocks()
    {
        return mAlignment;
    }

    public AlignedAssembly flipStrand()
    {
        final String assembly = SequenceUtil.reverseComplement(Assembly);
        final List<Alignment> flippedAlignments = new ArrayList<>();
        
        for(Alignment alignment : mAlignment)
        {
            final int startIndex = alignment.SequenceStartPosition - 1;
            final int newEndIndexExcl = assembly.length() - startIndex;
            final int newStartIndex = newEndIndexExcl - alignment.Length;

            if(alignment.isUnmapped())
            {
                flippedAlignments.add(new Alignment(alignment.Chromosome, 0,
                        newStartIndex + 1,
                        alignment.Length, false, 0));
            }
            else
            {
                flippedAlignments.add(new Alignment(
                        alignment.Chromosome, alignment.ReferenceStartPosition,
                        newStartIndex + 1,
                        alignment.Length, !alignment.Inverted, alignment.Quality));
            }
        }

        Collections.reverse(flippedAlignments);
        final AlignedAssembly flipped = new AlignedAssembly(assembly, Source.flipStrand(), flippedAlignments);

        for(ReadSupport support : readSupport())
        {
            int initialReadLength = support.Read.getLength();
            Read flippedRead = flipRead(support.Read);
            
            flipped.addEvidenceAt(flippedRead, getLength() - support.Index - initialReadLength);
        }

        flipped.recalculateBaseQuality();
        return flipped;
    }
}
