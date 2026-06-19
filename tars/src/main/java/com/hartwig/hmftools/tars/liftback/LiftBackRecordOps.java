package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.liftback.SaTagRewriter.SA_ATTRIBUTE;

import java.util.stream.Collectors;

import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.SequenceUtil;

// Pure per-record mutation primitives shared by the standalone SpliceLiftBack main and the REDUX
// LiftBackGroupProcessor. No group/orchestration state -- each method maps a resolved LiftBackResult onto
// a single SAMRecord (or derives the mate-cache view of it). Extracted from SpliceLiftBack so the
// concurrent REDUX path and the standalone path share one implementation.
public final class LiftBackRecordOps
{
    static final String XA_TAG = "XA";

    private LiftBackRecordOps() { }

    // Supplementaries whose own lift failed are mirrored onto their primary's lifted coords (keeping
    // the 0x800 flag) rather than marked unmapped -- htsjdk's validator rejects 0x4+0x800.
    public static void applyResultToRecord(
            final SAMRecord record, final LiftBackResult result, final LiftedMateInfoCache liftedMateInfoCache)
    {
        switch(result.category())
        {
            case UNMAPPED:
                return;

            case LIFT_FAILED:
                if(record.isSecondaryOrSupplementary() && mirrorOwnPrimaryOntoFailedSupp(record, liftedMateInfoCache))
                    return;

                record.setReadUnmappedFlag(true);
                record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
                record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
                record.setMappingQuality(0);
                record.setAttribute(XA_TAG, null);
                record.setAttribute(SA_ATTRIBUTE, null);
                return;

            default:
                final Cigar liftedCigar = TextCigarCodec.decode(result.finalCigar());
                record.setReferenceName(result.finalChrom());
                record.setAlignmentStart(result.finalPos());
                record.setCigar(liftedCigar);
                record.setReadNegativeStrandFlag(result.negativeStrand());
                record.setMappingQuality(result.updatedMapq());
                record.setAttribute(XA_TAG, buildLiftedXa(result));
                // XS:A:+/- on spliced records: downstream RNA tools (Isofox) rely on
                // transcript strand for stranded junction interpretation. bwa-mem2 emits XS:i: as
                // the sub-optimal alignment score on the same tag name -- clear first so the SAM
                // tag-type bookkeeping isn't ambiguous. Only set XS:A when the lifted cigar has
                // an N AND we know the transcript strand (i.e. the primary came off a tx contig).
                // Ref-only N-cigars from rescue/tail-extend don't have a strand threaded yet --
                // they ship without XS rather than risk a wrong call.
                record.setAttribute("XS", null);
                if(result.hasNCigar() && result.transcriptStrand() != 0)
                    record.setAttribute("XS", result.transcriptStrand() > 0 ? Character.valueOf('+') : Character.valueOf('-'));
        }
    }

    private static boolean mirrorOwnPrimaryOntoFailedSupp(
            final SAMRecord record, final LiftedMateInfoCache liftedMateInfoCache)
    {
        if(!record.getReadPairedFlag())
            return false;

        final LiftedMateInfo ownPrimary = liftedMateInfoCache.getOwnPrimary(record.getReadName(), record.getFirstOfPairFlag());
        if(ownPrimary == null || ownPrimary.unmapped())
            return false;

        record.setReferenceName(ownPrimary.chromosome());
        record.setAlignmentStart(ownPrimary.alignmentStart());
        record.setCigarString(ownPrimary.liftedCigar());
        record.setMappingQuality(0);
        record.setAttribute(XA_TAG, null);
        // SA stays -- rewriteSaTag translates tx-contig entries downstream and FragmentCoords requires it on supps
        return true;
    }

    public static LiftedMateInfo toLiftedMateInfo(final SAMRecord record, final LiftBackResult result)
    {
        if(willBeUnmapped(result))
            return LiftedMateInfo.UNMAPPED;

        final Cigar liftedCigar = TextCigarCodec.decode(result.finalCigar());
        final int alignmentEnd = result.finalPos() + liftedCigar.getReferenceLength() - 1;
        return LiftedMateInfo.mapped(result.finalChrom(), result.finalPos(), alignmentEnd, result.finalCigar(), result.negativeStrand());
    }

    // Mark a previously-mapped primary as unmapped, clearing every per-record tag the SAM spec invariant
    // would otherwise leave inconsistent. SA/XA/NH/MC all reference now-stale coords; ProperPair and
    // TLEN are meaningless once one end is unmapped.
    public static void markPrimaryUnmapped(final SAMRecord record)
    {
        record.setReadUnmappedFlag(true);
        record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
        record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
        record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
        record.setMappingQuality(0);
        record.setAttribute(SA_ATTRIBUTE, null);
        record.setAttribute(XA_TAG, null);
        record.setAttribute("NH", null);
        record.setAttribute("MC", null);
        if(record.getReadPairedFlag())
        {
            record.setProperPairFlag(false);
            record.setInferredInsertSize(0);
        }
    }

    // Single source of truth for whether a primary's final state ends up unmapped, used both at cache-
    // build time (so partner records see the correct MateUnmapped) and at apply time (where we mutate
    // the record itself).
    public static boolean willBeUnmapped(final LiftBackResult result)
    {
        if(result.category() == LiftBackCategory.UNMAPPED || result.category() == LiftBackCategory.LIFT_FAILED)
            return true;
        return exceedsMappingCap(result);
    }

    // bwa caps the XA alternative-hit list (run with -h 75 here); a read that maps to more loci than the cap
    // is emitted with MAPQ 0 and NO XA tag (the list is suppressed, not truncated). Such a read maps to too
    // many places to trust, so it is unmapped. inputMapq 0 distinguishes it from a unique read (MAPQ 60, no XA)
    // and numXaAlts 0 from an ordinary few-way multimapper (MAPQ 0 with its alts listed in XA).
    //
    // The cap is meaningful only for a GENOMIC (REF_ONLY) primary: there, 75+ suppressed hits are 75+ distinct
    // genomic loci. A tx-contig primary instead hits 75+ transcript contigs of one gene (a shared exon), which
    // all lift back to a single genomic locus, so its suppressed XA must NOT be read as "too many genomic
    // places". Gating on REF_ONLY keeps genuine genomic repeats unmapped while letting tx reads lift normally.
    public static boolean exceedsMappingCap(final LiftBackResult result)
    {
        return result.inputMapq() == 0 && result.numXaAlts() == 0 && result.comp() == LiftBackResult.Composition.REF_ONLY;
    }

    // The tx-contig MD/NM are stale once the read is lifted and (for rescue/extend/collapse/canon) recut.
    // MD is cigar-coupled and is dropped, not rebuilt: reconstructing it across a spliced read would need
    // reference spanning the whole intron. NM only needs the aligned M blocks (N/S/H contribute nothing),
    // so it is recomputed cheaply against the genomic reference. With no ref source (rescue + extend both
    // off) or on an unmapped record, NM is cleared as before.
    public static void refreshNmDropMd(final SAMRecord record, final RefSequenceSource refSource)
    {
        record.setAttribute("MD", null);

        if(record.getReadUnmappedFlag() || refSource == null)
        {
            record.setAttribute("NM", null);
            return;
        }

        final int editDistance = computeEditDistance(record, refSource);
        record.setAttribute("NM", editDistance >= 0 ? Integer.valueOf(editDistance) : null);
    }

    // NM per the SAM spec: mismatches in aligned blocks plus inserted and deleted bases. N (intron), soft
    // and hard clips, and padding contribute nothing, so each M block's reference is fetched on its own and
    // the intron is never read. Returns -1 if any block's reference is unavailable, so the caller leaves
    // NM unset rather than write a wrong value.
    private static int computeEditDistance(final SAMRecord record, final RefSequenceSource refSource)
    {
        final byte[] readBases = record.getReadBases();
        if(readBases == null || readBases.length == 0)
            return -1;

        final String chromosome = record.getReferenceName();
        int refPos = record.getAlignmentStart();
        int readIndex = 0;
        int editDistance = 0;

        for(final CigarElement element : record.getCigar().getCigarElements())
        {
            final int length = element.getLength();
            switch(element.getOperator())
            {
                case M:
                case EQ:
                case X:
                    final byte[] refBases = refSource.getBases(chromosome, refPos, refPos + length - 1);
                    if(refBases == null || refBases.length < length)
                        return -1;
                    for(int i = 0; i < length; ++i)
                    {
                        if(!SequenceUtil.basesEqual(readBases[readIndex + i], refBases[i]))
                            ++editDistance;
                    }
                    refPos += length;
                    readIndex += length;
                    break;

                case I:
                    editDistance += length;
                    readIndex += length;
                    break;

                case D:
                    editDistance += length;
                    refPos += length;
                    break;

                case N:
                    refPos += length;
                    break;

                case S:
                    readIndex += length;
                    break;

                case H:
                case P:
                    break;
            }
        }

        return editDistance;
    }

    public static String buildLiftedXa(final LiftBackResult result)
    {
        // Drop alts whose lifted span overlaps the primary choice's span: a tx-contig isoform read in a
        // shared exon lifts onto the primary's own coordinates (often a few bp off via softclip-boundary
        // or intron-rendering differences, so the exact-coord dedup in LiftBackResolver misses it). Those
        // entries re-report the locus bwa already placed the read at and carry no alternative-position
        // info. Anything non-overlapping (other locus, other chromosome, a different exon up the same gene)
        // is a genuine alternative mapping and is kept.
        final LiftedAlignment primary = result.liftedAlignments().stream()
                .filter(la -> la.IsPrimaryChoice)
                .findFirst().orElse(null);

        final String xa = result.liftedAlignments().stream()
                .filter(la -> !la.IsPrimaryChoice && !la.Dropped)
                .filter(la -> !overlapsPrimary(la, primary))
                .map(LiftBackRecordOps::formatXaEntry)
                .collect(Collectors.joining());
        return xa.isEmpty() ? null : xa;
    }

    private static boolean overlapsPrimary(final LiftedAlignment alt, final LiftedAlignment primary)
    {
        if(primary == null || !alt.LiftedChrom.equals(primary.LiftedChrom))
            return false;

        final int altStart = alt.LiftedPos;
        final int altEnd = altStart + TextCigarCodec.decode(alt.LiftedCigar).getReferenceLength() - 1;
        final int primStart = primary.LiftedPos;
        final int primEnd = primStart + TextCigarCodec.decode(primary.LiftedCigar).getReferenceLength() - 1;
        return altStart <= primEnd && altEnd >= primStart;
    }

    private static String formatXaEntry(final LiftedAlignment la)
    {
        return la.LiftedChrom + ',' + (la.ForwardStrand ? '+' : '-') + la.LiftedPos + ','
                + la.LiftedCigar + ',' + la.NumMismatches + ';';
    }
}
