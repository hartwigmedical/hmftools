package com.hartwig.hmftools.redux.splice;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public final class ContigTranslator
{
    public static TranslationResult translate(final ContigEntry contig, int contigPos, final Cigar contigCigar)
    {
        final List<BaseRegion> spans = contig.ExonSpans;
        if(spans.isEmpty() || contigPos < 1)
            return null;

        // locate the span that contains the read's leftmost mapped contig position
        int curSpanIdx = -1;
        int curGenomicPos = -1;
        int cumOffset = 0;
        for(int i = 0; i < spans.size(); ++i)
        {
            BaseRegion span = spans.get(i);
            if(contigPos <= cumOffset + span.baseLength())
            {
                curSpanIdx = i;
                curGenomicPos = span.start() + (contigPos - cumOffset - 1);
                break;
            }
            cumOffset += span.baseLength();
        }

        if(curSpanIdx < 0)
            return null; // read starts past the end of the contig

        int genomicStart = curGenomicPos;
        List<CigarElement> outElements = new ArrayList<>();
        List<BaseRegion> impliedIntrons = new ArrayList<>();

        for(CigarElement element : contigCigar.getCigarElements())
        {
            CigarOperator op = element.getOperator();
            int remaining = element.getLength();

            if(!op.consumesReferenceBases())
            {
                outElements.add(new CigarElement(remaining, op));
                continue;
            }

            while(remaining > 0)
            {
                BaseRegion currentSpan = spans.get(curSpanIdx);

                // if we've stepped past the current span, insert an intron and advance into the next span
                if(curGenomicPos > currentSpan.end())
                {
                    if(curSpanIdx + 1 >= spans.size())
                        return null; // read extends past last span — un-liftable

                    BaseRegion nextSpan = spans.get(curSpanIdx + 1);
                    int intronStart = currentSpan.end() + 1;
                    int intronEnd = nextSpan.start() - 1;

                    outElements.add(new CigarElement(intronEnd - intronStart + 1, CigarOperator.N));
                    impliedIntrons.add(new BaseRegion(intronStart, intronEnd));

                    ++curSpanIdx;
                    currentSpan = nextSpan;
                    curGenomicPos = currentSpan.start();
                }

                int remainInSpan = currentSpan.end() - curGenomicPos + 1;
                int take = Math.min(remaining, remainInSpan);

                outElements.add(new CigarElement(take, op));
                curGenomicPos += take;
                remaining -= take;
            }
        }

        return new TranslationResult(contig.Chromosome, genomicStart, new Cigar(outElements), impliedIntrons);
    }

    public static final class TranslationResult
    {
        public final String Chromosome;
        public final int GenomicStart;
        public final Cigar GenomicCigar;
        public final List<BaseRegion> ImpliedIntrons;

        public TranslationResult(
                final String chromosome, int genomicStart, final Cigar genomicCigar, final List<BaseRegion> impliedIntrons)
        {
            Chromosome = chromosome;
            GenomicStart = genomicStart;
            GenomicCigar = genomicCigar;
            ImpliedIntrons = impliedIntrons;
        }
    }
}
