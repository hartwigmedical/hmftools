package com.hartwig.hmftools.tars.liftback;

import java.util.HashSet;
import java.util.Set;

// Rewrites each SA tag entry from tx-contig to genomic coordinates.
// NM is not recomputed — it disagrees with the rewritten CIGAR once N operators are introduced, same as the read's own NM.
// Entries that fail to lift are dropped; survivors are deduped by (chrom, pos, strand, CIGAR).
public class SaTagRewriter
{
    public static final String SA_ATTRIBUTE = "SA";

    public static String rewriteSaTag(final String saTagValue, final LiftBackResolver resolver)
    {
        if(saTagValue == null || saTagValue.isEmpty())
            return null;

        final StringBuilder rewritten = new StringBuilder();
        final Set<String> seenEntryKeys = new HashSet<>();

        for(final String saEntry : saTagValue.split(";"))
        {
            if(saEntry.isEmpty())
                continue;

            final String[] parts = saEntry.split(",");
            if(parts.length < 6)
                continue;

            final String contig = parts[0];
            final int position;
            try
            {
                position = Integer.parseInt(parts[1]);
            }
            catch(NumberFormatException e)
            {
                continue;
            }
            final String strand = parts[2];
            final String cigarString = parts[3];
            final String mapQuality = parts[4];
            final String numMismatches = parts[5];

            final LiftedCoords lifted = resolver.liftCoords(contig, position, cigarString);
            if(lifted == null)
                continue;

            final String entryKey = lifted.chromosome() + ':' + lifted.position() + ':' + strand + ':' + lifted.cigarString();
            if(!seenEntryKeys.add(entryKey))
                continue;

            rewritten.append(lifted.chromosome()).append(',')
                    .append(lifted.position()).append(',')
                    .append(strand).append(',')
                    .append(lifted.cigarString()).append(',')
                    .append(mapQuality).append(',')
                    .append(numMismatches).append(';');
        }

        return rewritten.length() == 0 ? null : rewritten.toString();
    }
}
