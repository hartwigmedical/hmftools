package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.ALT_CONTIG_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;

// per-record lift-back driver. Two passes over the bwa-mem2 BAM aligned against ref + Tx contigs:
//  - pass 1: builds a LiftedMateInfoCache (one entry per paired primary record) so each read's
//            partner-mate fields can be rewritten in genomic coords.
//  - pass 2: resolves + applies + writes a coordinate-sorted BAM with mate fields patched.
// 1:1 record contract still holds: every input record produces exactly one output record + one TSV-A row + N TSV-B rows.
public class SpliceLiftBack
{
    private static final String XA_TAG = "XA";

    private final SpliceLiftBackConfig mConfig;
    private final ConfigBuilder mConfigBuilder;

    public SpliceLiftBack(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceLiftBackConfig(configBuilder);
        mConfigBuilder = configBuilder;
    }

    public void run()
    {
        final long startTimeMs = System.currentTimeMillis();

        final List<ContigEntry> contigEntries;
        if(mConfig.hasContigSidecar())
        {
            contigEntries = ContigSidecar.read(mConfig.ContigSidecarFile);
            RD_LOGGER.info("loaded {} contig entries", contigEntries.size());
        }
        else
        {
            contigEntries = List.of();
            RD_LOGGER.info("no contig sidecar supplied — running in pass-through mode (lift-back is a no-op)");
        }

        final LiftBackResolver resolver = new LiftBackResolver(contigEntries);
        validateBamAgainstSidecar(resolver.contigNames());

        final RrnaRegionIndex rrnaIndex = mConfig.FilterRrna
                ? RrnaRegionIndex.build(mConfigBuilder, RefGenomeVersion.from(mConfigBuilder))
                : null;
        final Set<String> filteredReadNames = collectFilteredReadNames(resolver, rrnaIndex);

        final LiftedMateInfoCache liftedMateInfoCache = buildLiftedMateInfoCache(resolver, filteredReadNames);

        emitLiftedBam(resolver, liftedMateInfoCache, filteredReadNames);

        RD_LOGGER.info("SpliceLiftBack complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // fail fast if the sidecar doesn't cover every _tx alt contig the BAM was aligned against.
    // a mismatch here was the root cause of an earlier silent-failure: alt contigs missing from the
    // segments map fell through as ref pass-through and leaked _tx names into mate fields downstream.
    private void validateBamAgainstSidecar(final Set<String> sidecarContigs)
    {
        Set<String> missing = new TreeSet<>();
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam)))
        {
            for(SAMSequenceRecord sq : reader.getFileHeader().getSequenceDictionary().getSequences())
            {
                String name = sq.getSequenceName();
                if(name.endsWith(ALT_CONTIG_SUFFIX) && !sidecarContigs.contains(name))
                    missing.add(name);
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read BAM header: " + mConfig.InputBam, e);
        }

        if(!missing.isEmpty())
        {
            throw new IllegalStateException(String.format(
                    "BAM/sidecar mismatch: %d alt contig(s) in BAM @SQ are absent from sidecar %s — first few: %s",
                    missing.size(), mConfig.ContigSidecarFile,
                    missing.stream().limit(5).collect(Collectors.joining(","))));
        }
    }

    // rRNA pre-pass: stream the BAM, resolve every record (primaries + supplementaries), flag any read
    // name where at least one record's post-lift coords overlap an rRNA gene region. Returned set drives
    // drop-by-read-pair in subsequent passes. Returns an empty set when filtering is disabled.
    private Set<String> collectFilteredReadNames(final LiftBackResolver resolver, final RrnaRegionIndex rrnaIndex)
    {
        if(rrnaIndex == null)
            return Set.of();

        final long startTimeMs = System.currentTimeMillis();
        final Set<String> filtered = new HashSet<>();
        int recordsScanned = 0;

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam)))
        {
            final SAMRecordIterator iter = samReader.iterator();
            while(iter.hasNext())
            {
                final SAMRecord record = iter.next();
                ++recordsScanned;

                if(record.getReadUnmappedFlag())
                    continue;

                final LiftBackResult result = resolver.resolve(record);
                if(result.Category == LiftBackCategory.UNMAPPED || result.Category == LiftBackCategory.LIFT_FAILED)
                    continue;

                final Cigar cigar = TextCigarCodec.decode(result.FinalCigar);
                final int liftedEnd = result.FinalPos + cigar.getReferenceLength() - 1;
                if(rrnaIndex.overlaps(result.FinalChrom, result.FinalPos, liftedEnd))
                    filtered.add(record.getReadName());
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("rRNA pre-pass BAM read failed", e);
            throw new RuntimeException(e);
        }

        RD_LOGGER.info("rRNA pre-pass scanned {} records, flagged {} read name(s), mins({})",
                recordsScanned, filtered.size(), runTimeMinsStr(startTimeMs));
        return filtered;
    }

    // pass 1: stream the input BAM, resolve each paired primary record, cache its lifted info keyed by read name
    // so pass 2 can patch the partner read's mate fields. Records belonging to rRNA-filtered read names are
    // skipped — they won't make the output BAM, and caching their mate info would waste memory.
    private LiftedMateInfoCache buildLiftedMateInfoCache(final LiftBackResolver resolver, final Set<String> filteredReadNames)
    {
        final long startTimeMs = System.currentTimeMillis();
        final LiftedMateInfoCache cache = new LiftedMateInfoCache();
        int primaryRecordsCached = 0;

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam)))
        {
            final SAMRecordIterator iter = samReader.iterator();
            while(iter.hasNext())
            {
                final SAMRecord record = iter.next();

                if(!record.getReadPairedFlag())
                    continue;
                if(record.isSecondaryOrSupplementary())
                    continue;
                if(filteredReadNames.contains(record.getReadName()))
                    continue;

                final LiftBackResult result = resolver.resolve(record);
                final LiftedMateInfo info = toLiftedMateInfo(record, result);
                cache.recordPrimaryAlignment(record.getReadName(), record.getFirstOfPairFlag(), info);
                ++primaryRecordsCached;
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("pass 1 BAM read failed", e);
            throw new RuntimeException(e);
        }

        RD_LOGGER.info("pass 1 cached {} primary alignments across {} read pairs, mins({})",
                primaryRecordsCached, cache.size(), runTimeMinsStr(startTimeMs));
        return cache;
    }

    // pass 2: stream the input BAM again, apply the lift-back to each record, patch mate fields from the cache,
    // and emit to a coordinate-sorted, indexed BAM. TSVs are written here so MateNum reflects firstOfPair.
    // Records whose read name is in filteredReadNames are dropped from the BAM (rRNA pair filter) but a TSV-A
    // audit row is still emitted with Filtered=true so the 1:1 record audit trail is preserved.
    private void emitLiftedBam(
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final Set<String> filteredReadNames)
    {
        final long startTimeMs = System.currentTimeMillis();
        final LiftBackStats stats = new LiftBackStats();
        int rrnaRecordsDropped = 0;

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.InputBam));
            LiftBackWriter writer = new LiftBackWriter(mConfig.formTsvAFile(), mConfig.formTsvBFile()))
        {
            final SAMFileHeader header = samReader.getFileHeader().clone();
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            final SAMFileWriterFactory factory = new SAMFileWriterFactory().setCreateIndex(true);

            try(SAMFileWriter samWriter = factory.makeBAMWriter(header, false, new File(mConfig.formOutputBam())))
            {
                final SAMRecordIterator iter = samReader.iterator();
                while(iter.hasNext())
                {
                    final SAMRecord record = iter.next();
                    final boolean dropForRrna = filteredReadNames.contains(record.getReadName());

                    final LiftBackResult result = resolver.resolve(record);
                    stats.record(record, result);

                    if(dropForRrna)
                    {
                        ++rrnaRecordsDropped;
                        writer.write(record, result, true);
                        continue;
                    }

                    applyResultToRecord(record, result, liftedMateInfoCache);
                    record.setAttribute(SA_ATTRIBUTE, rewriteSaTag(record.getStringAttribute(SA_ATTRIBUTE), resolver));
                    patchMateFields(record, liftedMateInfoCache);
                    samWriter.addAlignment(record);

                    writer.write(record, result, false);
                }
            }

            stats.writeSummary(mConfig.formSummaryFile());
        }
        catch(IOException e)
        {
            RD_LOGGER.error("pass 2 BAM/TSV I/O failed", e);
            throw new RuntimeException(e);
        }

        stats.logSummary();
        if(!filteredReadNames.isEmpty())
        {
            RD_LOGGER.info("rRNA filter dropped {} record(s) across {} read pair(s)",
                    rrnaRecordsDropped, filteredReadNames.size());
        }
        RD_LOGGER.info("pass 2 emit complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    // mutate the SAMRecord in place to reflect the LiftBackResult: lifted chrom/pos/CIGAR, or unmapped flag
    // for UNMAPPED/LIFT_FAILED. XA tag is rewritten from the deduped, lifted XA alts (Stage 4).
    //
    // Supplementary records whose own lift failed are mirrored onto their primary's lifted coordinates
    // (keeping the supplementary flag set) rather than being marked unmapped. htsjdk's BAM validator
    // rejects records with both the unmapped and supplementary flags set, and we never drop a read.
    static void applyResultToRecord(
            final SAMRecord record, final LiftBackResult result, final LiftedMateInfoCache liftedMateInfoCache)
    {
        switch(result.Category)
        {
            case UNMAPPED:
                // already unmapped; nothing to do
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
                Cigar liftedCigar = TextCigarCodec.decode(result.FinalCigar);
                record.setReferenceName(result.FinalChrom);
                record.setAlignmentStart(result.FinalPos);
                record.setCigar(liftedCigar);
                record.setReadNegativeStrandFlag(result.NegativeStrand);
                record.setMappingQuality(result.UpdatedMapq);
                record.setAttribute(XA_TAG, buildLiftedXa(result));
        }
    }

    // mirror a supplementary record onto its primary's lifted alignment. Returns false (caller falls back
    // to the unmap path) if the primary is missing from the cache or itself failed to lift.
    private static boolean mirrorOwnPrimaryOntoFailedSupp(
            final SAMRecord record, final LiftedMateInfoCache liftedMateInfoCache)
    {
        if(!record.getReadPairedFlag())
            return false;

        LiftedMateInfo ownPrimary = liftedMateInfoCache.getOwnPrimary(record.getReadName(), record.getFirstOfPairFlag());
        if(ownPrimary == null || ownPrimary.Unmapped)
            return false;

        record.setReferenceName(ownPrimary.Chromosome);
        record.setAlignmentStart(ownPrimary.AlignmentStart);
        record.setCigarString(ownPrimary.LiftedCigar);
        record.setMappingQuality(0);
        record.setAttribute(XA_TAG, null);
        // leave SA tag in place — rewriteSaTag translates any tx-contig entries to genomic coords downstream,
        // and redux's FragmentCoords requires SA on supplementary records.
        return true;
    }

    // derive the cached lifted mate info for a primary record from its resolved LiftBackResult.
    // UNMAPPED/LIFT_FAILED primaries are cached as unmapped so the partner can clear proper-pair on patch.
    static LiftedMateInfo toLiftedMateInfo(final SAMRecord record, final LiftBackResult result)
    {
        if(result.Category == LiftBackCategory.UNMAPPED || result.Category == LiftBackCategory.LIFT_FAILED)
            return LiftedMateInfo.unmapped();

        Cigar liftedCigar = TextCigarCodec.decode(result.FinalCigar);
        int alignmentEnd = result.FinalPos + liftedCigar.getReferenceLength() - 1;
        return LiftedMateInfo.mapped(result.FinalChrom, result.FinalPos, alignmentEnd, result.FinalCigar, result.NegativeStrand);
    }

    // build the bwa XA tag string from the lifted, deduped XA alts. Returns null when there are none, which
    // tells htsjdk to drop the tag entirely.
    static String buildLiftedXa(final LiftBackResult result)
    {
        // emit every alignment in the lifted set except (a) the chosen BAM primary, (b) any losers the
        // discriminator marked Dropped. Filtering by IsPrimaryChoice (rather than Source==SELF) is what
        // lets a swapped original-self appear in XA as an informative paralog hit.
        String xa = result.LiftedAlignments.stream()
                .filter(la -> !la.IsPrimaryChoice && !la.Dropped)
                .map(SpliceLiftBack::formatXaEntry)
                .collect(Collectors.joining());
        return xa.isEmpty() ? null : xa;
    }

    private static String formatXaEntry(final LiftedAlignment la)
    {
        return la.LiftedChrom + ',' + (la.ForwardStrand ? '+' : '-') + la.LiftedPos + ','
                + la.LiftedCigar + ',' + la.NumMismatches + ';';
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceLiftBackConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SpliceLiftBack liftBack = new SpliceLiftBack(configBuilder);
        liftBack.run();
    }
}
