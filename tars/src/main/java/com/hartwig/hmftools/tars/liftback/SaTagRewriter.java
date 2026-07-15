package com.hartwig.hmftools.tars.liftback;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

// Rewrites each SA tag entry from tx-contig to genomic coordinates.
// NM is not recomputed - it disagrees with the rewritten CIGAR once N operators are introduced, same as the read's own NM.
// Entries that fail to lift are dropped; survivors are deduped by (chrom, pos, strand, CIGAR).
public class SaTagRewriter
{
    public static String rewriteSaTag(final String saTagValue, final LiftBackResolver resolver)
    {
        return rewriteSaTag(saTagValue, resolver, Collections.emptySet());
    }

    // excludeKeys: lifted entry keys (chrom:pos:strand:cigar, clips normalised H->S) of supplementaries dropped from
    // the output (excluded region / orphan / low-AS). Their SA entries are removed so the primary never references a
    // supp that isn't emitted. Clips are normalised because the dropped supp RECORD is hard-clipped while its SA-tag
    // entry here is soft-clipped - see LiftBackGroupProcessor.suppSaKey.
    public static String rewriteSaTag(final String saTagValue, final LiftBackResolver resolver, final Set<String> excludeKeys)
    {
        if(saTagValue == null || saTagValue.isEmpty())
        {
            return null;
        }

        List<SupplementaryReadData> alignments = SupplementaryReadData.extractAlignments(saTagValue);
        if(alignments == null)
        {
            return null;
        }

        StringBuilder rewritten = new StringBuilder();
        Set<String> seenEntryKeys = new HashSet<>();

        for(SupplementaryReadData alignment : alignments)
        {
            LiftedCoords lifted = resolver.liftCoords(alignment.Chromosome, alignment.Position, alignment.Cigar);
            if(lifted == null)
                continue;

            String entryKey = liftedEntryKey(lifted, alignment.Strand);
            if(excludeKeys.contains(entryKey))
                continue;
            if(!seenEntryKeys.add(entryKey))
                continue;

            SupplementaryReadData liftedData = new SupplementaryReadData(
                    lifted.chromosome(), lifted.position(), alignment.Strand, lifted.cigarString(),
                    alignment.MapQuality, alignment.NM);
            rewritten.append(liftedData.asSamTag()).append(';');
        }

        return rewritten.length() == 0 ? null : rewritten.toString();
    }

    // Lifted SA-entry key (chrom:pos:strand:cigar, clips normalised H->S). Shared by rewriteSaTag and the
    // dropped-supp exclude set (LiftBackGroupProcessor) so both keys are built from the SAME lifted coords and
    // can never drift apart - a dropped supp's own lift path (micro-anchor trimmed) would otherwise produce a
    // different cigar/pos than lifting its SA entry here, leaving the primary's SA referencing a missing supp.
    static String liftedEntryKey(final LiftedCoords lifted, final char strand)
    {
        return lifted.chromosome() + ':' + lifted.position() + ':' + strand + ':' + lifted.cigarString().replace('H', 'S');
    }
}
