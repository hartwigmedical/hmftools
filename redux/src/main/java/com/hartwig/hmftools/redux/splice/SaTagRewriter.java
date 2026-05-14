package com.hartwig.hmftools.redux.splice;

import java.util.HashSet;
import java.util.Set;

// rewrites the SA tag of a SAMRecord so each chimeric alignment entry (primary <-> supplementaries)
// references genomic coordinates rather than transcript-contig coordinates.
//
// SA tag entry format: rname,pos,strand,CIGAR,mapQ,NM; (semicolon-separated, htsjdk SAMTag.SA).
// Strand char and mapQ pass through unchanged. NM is intentionally not recomputed (would otherwise
// disagree with the rewritten CIGAR when junction N operators are introduced) — same limitation as
// the read's own NM tag in SpliceLiftBack.
//
// Entries that fail to lift (alt-contig position outside any transcript span, malformed entry) are
// dropped. Remaining entries are deduped by lifted (chromosome, position, strand, CIGAR).
public class SaTagRewriter
{
    public static final String SA_ATTRIBUTE = "SA";

    public static String rewriteSaTag(final String saTagValue, final LiftBackResolver resolver)
    {
        if(saTagValue == null || saTagValue.isEmpty())
            return null;

        StringBuilder rewritten = new StringBuilder();
        Set<String> seenEntryKeys = new HashSet<>();

        for(String saEntry : saTagValue.split(";"))
        {
            if(saEntry.isEmpty())
                continue;

            String[] saEntryParts = saEntry.split(",");
            if(saEntryParts.length < 6)
                continue;

            String contig = saEntryParts[0];
            int position;
            try
            {
                position = Integer.parseInt(saEntryParts[1]);
            }
            catch(NumberFormatException e)
            {
                continue;
            }
            String strand = saEntryParts[2];
            String cigarString = saEntryParts[3];
            String mapQuality = saEntryParts[4];
            String numMismatches = saEntryParts[5];

            LiftBackResolver.LiftedCoords liftedCoords = resolver.liftCoords(contig, position, cigarString);
            if(liftedCoords == null)
                continue;

            String entryKey = liftedCoords.Chromosome + ':' + liftedCoords.Position + ':' + strand
                    + ':' + liftedCoords.CigarString;
            if(!seenEntryKeys.add(entryKey))
                continue;

            rewritten.append(liftedCoords.Chromosome).append(',')
                    .append(liftedCoords.Position).append(',')
                    .append(strand).append(',')
                    .append(liftedCoords.CigarString).append(',')
                    .append(mapQuality).append(',')
                    .append(numMismatches).append(';');
        }

        return rewritten.length() == 0 ? null : rewritten.toString();
    }
}
