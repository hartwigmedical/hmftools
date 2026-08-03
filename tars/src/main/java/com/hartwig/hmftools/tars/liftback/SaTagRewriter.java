package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_NEG_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.SAMRecord;

// Rewrites each SA tag entry from tx-contig to genomic coordinates. NM is left as bwa reported it: once an N operator
// is introduced it no longer agrees with the rewritten CIGAR, exactly as for the read's own NM. Entries that fail to
// lift are dropped, and survivors are deduped on the lifted key.
public class SaTagRewriter
{
    public static String rewriteSaTag(final String saTagValue, final ContigTranslator translator)
    {
        return rewriteSaTag(saTagValue, translator, Collections.emptySet());
    }

    // excludeKeys holds the lifted keys of supplementaries that will not be emitted (excluded region, orphan, low
    // score), so a primary is never left pointing at a record that is absent from the BAM.
    public static String rewriteSaTag(
            final String saTagValue, final ContigTranslator translator, final Set<String> excludeKeys)
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
            LiftedAlignment lifted = translator.liftAlignment(
                    alignment.Chromosome, alignment.Position, alignment.Cigar, alignment.NM, alignment.isForwardOrient());

            if(lifted == null)
            {
                continue;
            }

            String entryKey = liftedEntryKey(lifted);
            if(excludeKeys.contains(entryKey) || !seenEntryKeys.add(entryKey))
            {
                continue;
            }

            SupplementaryReadData liftedData = new SupplementaryReadData(
                    lifted.LiftedChromosome, lifted.LiftedPos, alignment.Strand, lifted.LiftedCigar,
                    alignment.MapQuality, alignment.NM);

            rewritten.append(liftedData.asSamTag()).append(';');
        }

        return rewritten.length() == 0 ? null : rewritten.toString();
    }

    // Clips are normalised H->S because a dropped supplementary's own record is hard-clipped while the SA entry
    // describing it is soft-clipped, and both have to key the same way.
    private static String liftedEntryKey(final LiftedAlignment lifted)
    {
        return lifted.LiftedChromosome + ':' + lifted.LiftedPos
                + ':' + (lifted.ForwardStrand ? SUPP_POS_STRAND : SUPP_NEG_STRAND)
                + ':' + lifted.LiftedCigar.replace('H', 'S');
    }

    // Lifted keys of the supplementaries that will not be emitted, keyed the same way rewriteSaTag keys them: from the
    // supplementary's own bwa placement, which the record still holds at this point. Keying off the post-liftback
    // result instead would use a trimmed CIGAR that never matches, leaving the primary's SA dangling.
    static Set<String> droppedSuppSaKeys(
            final List<SAMRecord> records, final boolean[] willEmit, final SAMRecord primary,
            final ContigTranslator translator)
    {
        Set<String> keys = new HashSet<>();

        for(int i = 0; i < records.size(); ++i)
        {
            SAMRecord record = records.get(i);
            if(willEmit[i] || record == primary || !record.getSupplementaryAlignmentFlag())
            {
                continue;
            }

            LiftedAlignment lifted = translator.liftAlignment(
                    record.getReferenceName(), record.getAlignmentStart(), record.getCigarString(),
                    0, !record.getReadNegativeStrandFlag());

            if(lifted != null)
            {
                keys.add(liftedEntryKey(lifted));
            }
        }

        return keys;
    }

    // When the discriminator moved the primary off bwa's placement, a surviving supplementary's own SA still points at
    // bwa's pre-move locus, which is never emitted. This rebuilds it to reference the emitted primary followed by any
    // sibling emitted supplementaries, so downstream FragmentCoords reads the right primary coords.
    static String buildSwappedSuppSa(
            final List<SAMRecord> records, final LiftedRecord[] resolved, final boolean[] willEmit,
            final LiftedRecord primaryResult, final int suppIndex)
    {
        if(!primaryResult.hasPlacement())
        {
            return null;
        }

        StringBuilder saTag = new StringBuilder(saEntry(primaryResult));

        for(int i = 0; i < records.size(); ++i)
        {
            if(i == suppIndex || !willEmit[i] || !records.get(i).getSupplementaryAlignmentFlag())
            {
                continue;
            }

            LiftedRecord suppResult = resolved[i];
            if(suppResult == null || !suppResult.hasPlacement()
                    || suppResult.finalCigar().equals(SAMRecord.NO_ALIGNMENT_CIGAR))
            {
                continue;
            }

            saTag.append(saEntry(suppResult));
        }

        return saTag.length() > 0 ? saTag.toString() : null;
    }

    // one SA entry for a lifted placement; NM is left 0 since downstream keys on the coordinates
    private static String saEntry(final LiftedRecord result)
    {
        return result.finalChromosome() + ',' + result.finalPos()
                + ',' + (result.negativeStrand() ? SUPP_NEG_STRAND : SUPP_POS_STRAND)
                + ',' + result.finalCigar() + ',' + result.updatedMapQuality() + ",0;";
    }
}
