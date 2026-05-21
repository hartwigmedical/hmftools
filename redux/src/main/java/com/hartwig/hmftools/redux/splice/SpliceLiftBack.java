package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.redux.ReduxConfig.APP_NAME;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.MateFieldPatcher.patchMateFields;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.SA_ATTRIBUTE;
import static com.hartwig.hmftools.redux.splice.SaTagRewriter.rewriteSaTag;
import static com.hartwig.hmftools.redux.splice.SpliceCommon.ALT_CONTIG_SUFFIX;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;

// Two passes over the bwa-mem2 BAM aligned against ref + tx contigs:
//   pass 1: build a LiftedMateInfoCache so each record's partner-mate fields can be rewritten in genomic coords.
//   pass 2: resolve + apply + write a coordinate-sorted BAM with mate fields patched.
// 1:1 record contract: every input record produces exactly one output record + one TSV-A row + N TSV-B rows.
public class SpliceLiftBack
{
    private static final String XA_TAG = "XA";

    private final SpliceLiftBackConfig mConfig;

    // in -emit_secondaries mode this is a name-sorted intermediate so primary + 0x100 secondaries are
    // contiguous in the stream. Otherwise it equals mConfig.InputBam.
    private String mWorkingInputBam;

    public SpliceLiftBack(final ConfigBuilder configBuilder)
    {
        mConfig = new SpliceLiftBackConfig(configBuilder);
    }

    public void run()
    {
        final long startTimeMs = System.currentTimeMillis();

        mWorkingInputBam = mConfig.EmitSecondaries ? nameSortInput() : mConfig.InputBam;

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

        final LiftBackResolver resolver = new LiftBackResolver(contigEntries, mConfig.EmitSecondaries);
        validateBamAgainstSidecar(resolver.contigNames());

        final LiftedMateInfoCache liftedMateInfoCache = buildLiftedMateInfoCache(resolver);

        emitLiftedBam(resolver, liftedMateInfoCache);

        if(mConfig.EmitSecondaries && !mWorkingInputBam.equals(mConfig.InputBam))
        {
            try
            {
                Files.deleteIfExists(Paths.get(mWorkingInputBam));
            }
            catch(IOException e)
            {
                RD_LOGGER.warn("could not delete name-sorted intermediate {}: {}", mWorkingInputBam, e.toString());
            }
        }

        RD_LOGGER.info("SpliceLiftBack complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private String nameSortInput()
    {
        if(mConfig.BamToolPath == null)
            throw new IllegalStateException("-emit_secondaries requires -bamtool for the name-sort pre-pass");

        final String nameSortedBam = mConfig.OutputDir + "splice_lift_back.name_sorted_input.bam";
        RD_LOGGER.info("name-sorting input via {}: {}", fromPath(mConfig.BamToolPath), nameSortedBam);

        final long startTimeMs = System.currentTimeMillis();
        final java.util.List<String> command = new ArrayList<>();
        command.add(mConfig.BamToolPath);
        command.add("sort");
        command.add("-n");
        if(mConfig.Threads > 1)
        {
            command.add("-@");
            command.add(String.valueOf(mConfig.Threads));
        }
        if(fromPath(mConfig.BamToolPath) == BamToolName.SAMTOOLS)
        {
            command.add("-O");
            command.add("bam");
        }
        command.add(mConfig.InputBam);
        command.add("-o");
        command.add(nameSortedBam);

        try
        {
            final Process process = new ProcessBuilder(command).redirectErrorStream(true).start();
            final int exitCode = process.waitFor();
            if(exitCode != 0)
                throw new RuntimeException("name-sort exit code " + exitCode);
        }
        catch(InterruptedException | IOException e)
        {
            throw new RuntimeException("name-sort failed: " + e, e);
        }

        RD_LOGGER.info("name-sort complete, mins({})", runTimeMinsStr(startTimeMs));
        return nameSortedBam;
    }

    // fails fast on a BAM/sidecar mismatch — root cause of a prior silent failure where alt contigs
    // missing from the sidecar fell through as ref pass-through and leaked _tx names into mate fields.
    private void validateBamAgainstSidecar(final Set<String> sidecarContigs)
    {
        Set<String> missing = new TreeSet<>();
        try(SamReader reader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mWorkingInputBam)))
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
            throw new RuntimeException("failed to read BAM header: " + mWorkingInputBam, e);
        }

        if(!missing.isEmpty())
        {
            throw new IllegalStateException(String.format(
                    "BAM/sidecar mismatch: %d alt contig(s) in BAM @SQ are absent from sidecar %s — first few: %s",
                    missing.size(), mConfig.ContigSidecarFile,
                    missing.stream().limit(5).collect(Collectors.joining(","))));
        }
    }

    private LiftedMateInfoCache buildLiftedMateInfoCache(final LiftBackResolver resolver)
    {
        final long startTimeMs = System.currentTimeMillis();
        final LiftedMateInfoCache cache = new LiftedMateInfoCache();
        int primaryRecordsCached = 0;

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mWorkingInputBam)))
        {
            final SAMRecordIterator iter = samReader.iterator();
            while(iter.hasNext())
            {
                final SAMRecord record = iter.next();

                if(!record.getReadPairedFlag())
                    continue;
                if(record.isSecondaryOrSupplementary())
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

    private void emitLiftedBam(final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache)
    {
        final long startTimeMs = System.currentTimeMillis();
        final LiftBackStats stats = new LiftBackStats();

        final String unsortedBam = mConfig.formUnsortedBam();

        try(SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mWorkingInputBam));
            LiftBackWriter writer = new LiftBackWriter(mConfig.formTsvAFile(), mConfig.formTsvBFile()))
        {
            final SAMFileHeader header = samReader.getFileHeader().clone();
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            if(mConfig.hasContigSidecar())
                stripAltContigsFromHeader(header);

            try(SAMFileWriter samWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(unsortedBam)))
            {
                final SAMRecordIterator iter = samReader.iterator();

                if(mConfig.EmitSecondaries)
                    streamGrouped(iter, samWriter, writer, resolver, liftedMateInfoCache, stats);
                else
                    streamPerRecord(iter, samWriter, writer, resolver, liftedMateInfoCache, stats);
            }

            stats.writeSummary(mConfig.formSummaryFile());
        }
        catch(IOException e)
        {
            RD_LOGGER.error("pass 2 BAM/TSV I/O failed", e);
            throw new RuntimeException(e);
        }

        stats.logSummary();
        RD_LOGGER.info("pass 2 emit complete, mins({})", runTimeMinsStr(startTimeMs));

        sortAndIndex(unsortedBam, mConfig.formOutputBam());
    }

    private void streamPerRecord(
            final SAMRecordIterator iter, final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        while(iter.hasNext())
        {
            final SAMRecord record = iter.next();
            final LiftBackResult result = resolver.resolve(record);
            final int nh = countKeptForNh(result);
            applyAndWriteRecord(record, result, nh, samWriter, writer, resolver, liftedMateInfoCache, stats);
        }
    }

    // -emit_secondaries mode on a name-sorted input: buffer records sharing a read name and fold
    // secondaries into the primary's alignment set so the discriminator sees the full alt picture.
    // Supplementaries stay independent — they're split-read components, not alts.
    private void streamGrouped(
            final SAMRecordIterator iter, final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        final List<SAMRecord> buffer = new ArrayList<>();
        String currentName = null;
        while(iter.hasNext())
        {
            final SAMRecord record = iter.next();
            final String name = record.getReadName();
            if(currentName != null && !name.equals(currentName))
            {
                flushReadNameGroup(buffer, samWriter, writer, resolver, liftedMateInfoCache, stats);
                buffer.clear();
            }
            buffer.add(record);
            currentName = name;
        }
        if(!buffer.isEmpty())
        {
            flushReadNameGroup(buffer, samWriter, writer, resolver, liftedMateInfoCache, stats);
        }
    }

    private void flushReadNameGroup(
            final List<SAMRecord> group, final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        final List<SAMRecord> firstOfPair = new ArrayList<>();
        final List<SAMRecord> secondOfPair = new ArrayList<>();
        for(final SAMRecord record : group)
        {
            if(!record.getReadPairedFlag() || record.getFirstOfPairFlag())
                firstOfPair.add(record);
            else
                secondOfPair.add(record);
        }

        emitMateGroup(firstOfPair, samWriter, writer, resolver, liftedMateInfoCache, stats);
        emitMateGroup(secondOfPair, samWriter, writer, resolver, liftedMateInfoCache, stats);
    }

    private void emitMateGroup(
            final List<SAMRecord> records, final SAMFileWriter samWriter, final LiftBackWriter writer,
            final LiftBackResolver resolver, final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        if(records.isEmpty())
            return;

        SAMRecord primary = null;
        final List<SAMRecord> secondaries = new ArrayList<>();
        for(final SAMRecord record : records)
        {
            if(record.getSupplementaryAlignmentFlag())
                continue;
            if(record.isSecondaryAlignment())
                secondaries.add(record);
            else if(primary == null)
                primary = record;
        }

        final boolean ranGroupResolve = primary != null && !primary.getReadUnmappedFlag();
        final LiftBackResult primaryResult;
        if(ranGroupResolve)
            primaryResult = resolver.resolvePrimaryWithSecondaries(primary, secondaries);
        else
            primaryResult = primary != null ? resolver.resolve(primary) : null;

        final LiftBackResult[] resolved = new LiftBackResult[records.size()];
        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            resolved[i] = (record == primary && primaryResult != null) ? primaryResult : resolver.resolve(record);
        }

        final Set<String> keptKeys = (mConfig.DropLosingAlts && ranGroupResolve)
                ? buildKeptKeys(primaryResult)
                : null;

        final boolean[] willEmit = new boolean[records.size()];
        int emittedNonSupp = 0;
        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            final LiftBackResult result = resolved[i];
            boolean drop = false;
            if(keptKeys != null && record != primary && record.isSecondaryAlignment())
            {
                final String key = liftedKey(result.finalChrom(), result.finalPos(), result.finalCigar());
                if(!keptKeys.contains(key))
                    drop = true;
            }
            willEmit[i] = !drop;
            if(willEmit[i] && !record.getSupplementaryAlignmentFlag())
                ++emittedNonSupp;
        }
        final int nh = Math.max(emittedNonSupp, 1);

        for(int i = 0; i < records.size(); i++)
        {
            final SAMRecord record = records.get(i);
            final LiftBackResult result = resolved[i];
            if(!willEmit[i])
            {
                stats.record(record, result);
                continue;
            }
            applyAndWriteRecord(record, result, nh, samWriter, writer, resolver, liftedMateInfoCache, stats);
        }
    }

    private static Set<String> buildKeptKeys(final LiftBackResult result)
    {
        final Set<String> keys = new HashSet<>();
        for(final LiftedAlignment alignment : result.liftedAlignments())
        {
            if(!alignment.Dropped)
                keys.add(liftedKey(alignment.LiftedChrom, alignment.LiftedPos, alignment.LiftedCigar));
        }
        return keys;
    }

    private static String liftedKey(final String chrom, final int pos, final String cigar)
    {
        return chrom + ":" + pos + ":" + cigar;
    }

    private static int countKeptForNh(final LiftBackResult result)
    {
        if(result.liftedAlignments().isEmpty())
            return 1;
        int count = 0;
        for(final LiftedAlignment alignment : result.liftedAlignments())
        {
            if(!alignment.Dropped)
                ++count;
        }
        return Math.max(count, 1);
    }

    private void applyAndWriteRecord(
            final SAMRecord record, final LiftBackResult result, final int nh, final SAMFileWriter samWriter,
            final LiftBackWriter writer, final LiftBackResolver resolver,
            final LiftedMateInfoCache liftedMateInfoCache, final LiftBackStats stats)
    {
        stats.record(record, result);

        applyResultToRecord(record, result, liftedMateInfoCache);
        record.setAttribute(SA_ATTRIBUTE, rewriteSaTag(record.getStringAttribute(SA_ATTRIBUTE), resolver));
        patchMateFields(record, liftedMateInfoCache);

        if(result.category() != LiftBackCategory.UNMAPPED)
            record.setAttribute("NH", nh);

        if(mConfig.StarMapqLadder && !record.getReadUnmappedFlag())
        {
            final int ladderMapq = NhMapqLadder.starLadder(nh);
            final int newMapq = result.category().txWon() ? ladderMapq : Math.min(ladderMapq, result.inputMapq());
            record.setMappingQuality(newMapq);
        }

        if(mConfig.UnmapAboveNh > 0 && nh > mConfig.UnmapAboveNh
                && !record.getReadUnmappedFlag()
                && !record.isSecondaryAlignment()
                && !record.getSupplementaryAlignmentFlag())
        {
            record.setReadUnmappedFlag(true);
            record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            record.setMappingQuality(0);
            record.setAttribute(SA_ATTRIBUTE, null);
            if(record.getReadPairedFlag())
            {
                record.setProperPairFlag(false);
                record.setInferredInsertSize(0);
            }
        }

        if(mConfig.UnmapBelowMapq > 0
                && !record.getReadUnmappedFlag()
                && !record.isSecondaryAlignment()
                && !record.getSupplementaryAlignmentFlag()
                && record.getMappingQuality() < mConfig.UnmapBelowMapq)
        {
            record.setReadUnmappedFlag(true);
            record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
            record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            record.setMappingQuality(0);
            record.setAttribute(SA_ATTRIBUTE, null);
            if(record.getReadPairedFlag())
            {
                record.setProperPairFlag(false);
                record.setInferredInsertSize(0);
            }
        }

        stripMdNm(record);
        samWriter.addAlignment(record);

        try
        {
            writer.write(record, result);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
    }

    // MD/NM from the original tx-contig alignment are stale post-lift; downstream tools recompute them.
    private static void stripMdNm(final SAMRecord record)
    {
        record.setAttribute("MD", null);
        record.setAttribute("NM", null);
    }

    private void sortAndIndex(final String unsortedBam, final String sortedBam)
    {
        if(mConfig.BamToolPath == null)
        {
            RD_LOGGER.info("no -{} configured; leaving unsorted BAM at {}", BamToolName.BAMTOOL_PATH, unsortedBam);
            return;
        }

        final BamToolName toolName = fromPath(mConfig.BamToolPath);
        RD_LOGGER.info("sorting BAM via {}: {}", toolName, sortedBam);

        if(!BamOperations.sortBam(toolName, mConfig.BamToolPath, unsortedBam, sortedBam, mConfig.Threads))
            throw new RuntimeException("failed to sort BAM: " + unsortedBam);

        if(toolName == BamToolName.SAMTOOLS)
        {
            if(!BamOperations.indexBam(toolName, mConfig.BamToolPath, sortedBam, mConfig.Threads))
                throw new RuntimeException("failed to index BAM: " + sortedBam);
        }

        try
        {
            Files.deleteIfExists(Paths.get(unsortedBam));
        }
        catch(IOException e)
        {
            RD_LOGGER.warn("could not delete intermediate {}: {}", unsortedBam, e.toString());
        }
    }

    private static void stripAltContigsFromHeader(final SAMFileHeader header)
    {
        final SAMSequenceDictionary dict = header.getSequenceDictionary();
        final List<SAMSequenceRecord> kept = new ArrayList<>();
        int dropped = 0;

        for(final SAMSequenceRecord seq : dict.getSequences())
        {
            if(seq.getSequenceName().endsWith(ALT_CONTIG_SUFFIX))
                ++dropped;
            else
                kept.add(seq);
        }

        header.setSequenceDictionary(new SAMSequenceDictionary(kept));
        RD_LOGGER.info("dropped {} alt contig @SQ entries from output header", dropped);
    }

    // Supplementaries whose own lift failed are mirrored onto their primary's lifted coords (keeping
    // the 0x800 flag) rather than marked unmapped — htsjdk's validator rejects 0x4+0x800.
    static void applyResultToRecord(
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
        // SA stays — rewriteSaTag translates tx-contig entries downstream and FragmentCoords requires it on supps
        return true;
    }

    static LiftedMateInfo toLiftedMateInfo(final SAMRecord record, final LiftBackResult result)
    {
        if(result.category() == LiftBackCategory.UNMAPPED || result.category() == LiftBackCategory.LIFT_FAILED)
            return LiftedMateInfo.UNMAPPED;

        final Cigar liftedCigar = TextCigarCodec.decode(result.finalCigar());
        final int alignmentEnd = result.finalPos() + liftedCigar.getReferenceLength() - 1;
        return LiftedMateInfo.mapped(result.finalChrom(), result.finalPos(), alignmentEnd, result.finalCigar(), result.negativeStrand());
    }

    static String buildLiftedXa(final LiftBackResult result)
    {
        final String xa = result.liftedAlignments().stream()
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
        final ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        SpliceLiftBackConfig.addConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        new SpliceLiftBack(configBuilder).run();
    }
}
