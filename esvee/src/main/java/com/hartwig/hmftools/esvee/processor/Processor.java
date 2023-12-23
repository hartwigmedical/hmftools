package com.hartwig.hmftools.esvee.processor;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.TimeUnit;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.Junction;
import com.hartwig.hmftools.esvee.RegionOfInterest;
import com.hartwig.hmftools.esvee.assembly.AssemblyExtender;
import com.hartwig.hmftools.esvee.assembly.JunctionMetrics;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssembler;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.Alignment;
import com.hartwig.hmftools.esvee.models.ExtendedAssembly;
import com.hartwig.hmftools.esvee.models.GappedAssembly;
import com.hartwig.hmftools.esvee.models.PrimaryAssembly;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.models.Sequence;
import com.hartwig.hmftools.esvee.models.SupportedAssembly;
import com.hartwig.hmftools.esvee.output.VCFWriter;
import com.hartwig.hmftools.esvee.output.VariantLine;
import com.hartwig.hmftools.esvee.output.html.SummaryPageGenerator;
import com.hartwig.hmftools.esvee.output.html.VariantCallPageGenerator;
import com.hartwig.hmftools.esvee.util.CSVReader;
import com.hartwig.hmftools.esvee.util.CSVWriter;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;
import com.hartwig.hmftools.esvee.util.ParallelMapper;
import com.hartwig.hmftools.esvee.util.RangeUtils;
import com.hartwig.hmftools.esvee.util.StringUtils;
import com.hartwig.hmftools.esvee.util.Timeout;
import com.hartwig.hmftools.esvee.ImmutableJunction;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.BAMStreamWriter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;

public class Processor
{
    private static final Logger LOGGER = LogManager.getLogger(Processor.class);

    private final Context mContext;
    private final HomologySlider mHomologySlider;

    public final OverallCounters Counters = new OverallCounters();

    public Processor(final Context context)
    {
        mContext = context;
        mHomologySlider = new HomologySlider(mContext.ReferenceGenome);
    }

    public List<VariantCall> run()
    {
        try
        {
            final List<ImmutableJunction> junctions =
                    new CSVReader<>(ImmutableJunction.class, mContext.Config.junctionFile().getAbsolutePath()).readToEnd();
            return run(junctions);
        }
        catch(final IOException exception)
        {
            throw new RuntimeException(exception);
        }
    }

    public List<VariantCall> run(final List<? extends Junction> junctions)
    {
        final long startTimeNanos = System.nanoTime();
        Counters.JunctionsProcessed.add(junctions.size());
        LOGGER.info("Starting primary assembly on {} junctions", junctions.size());

        // Primary Junction Assembly
        final List<PrimaryAssemblyResult> primaryAssemblyResults = ParallelMapper.mapWithProgress(
                "Primary Assembly", Counters.PrimaryAssemblyTime, mContext.Executor, junctions,
                junction -> PrimaryAssembler.process(mContext, junction, Counters.PrimaryAssemblerCounters));
        final int primaryAssemblyCount = primaryAssemblyResults.stream().mapToInt(r -> r.Assemblies.size()).sum();
        LOGGER.info("Created {} primary assemblies in {}", primaryAssemblyCount, Counters.PrimaryAssemblyTime.formatValue());
        if(!mContext.Config.debug())
            junctions.clear();

        // Inter-junction deduplication
        final List<PrimaryAssembly> primaryAssemblies = Counters.InterJunctionDeduplicationTime.time(
                () -> consolidateNearbyAssemblies(flatten(primaryAssemblyResults)));
        LOGGER.info("Reduced to {} assemblies in {}", primaryAssemblies.size(), Counters.InterJunctionDeduplicationTime.formatValue());
        if(!mContext.Config.debug())
            primaryAssemblyResults.clear();

        if(mContext.Config.dropGermline())
            primaryAssemblies.removeIf(assembly -> assembly.getSupportRecords().stream().anyMatch(Record::isGermline));

        // Assembly Extension
        final List<ExtendedAssembly> extendedAssemblies = ParallelMapper.mapWithProgress(
                        "Assembly Extension", Counters.ExtensionTime, mContext.Executor, primaryAssemblies,
                        assembly -> AssemblyExtender.process(mContext, assembly, Counters.AssemblyExtenderCounters)).stream()
                .flatMap(Collection::stream)
                .collect(Collectors.toList());
        if(!mContext.Config.debug())
            primaryAssemblies.clear();
        LOGGER.info("Created {} extended assemblies in {}", extendedAssemblies.size(), Counters.ExtensionTime.formatValue());

        // Primary phasing
        final List<Set<ExtendedAssembly>> primaryPhaseSets = Counters.PrimaryPhasingTime.time(
                () -> PrimaryPhasing.run(extendedAssemblies));
        LOGGER.info("Created {} primary phase sets in {}", primaryPhaseSets.size(), Counters.PrimaryPhasingTime.formatValue());
        if(!mContext.Config.debug())
            extendedAssemblies.clear();

        // Phased assembly merging
        final List<Set<ExtendedAssembly>> mergedPhaseSets = ParallelMapper.mapWithProgress(
                "Phased Merging", Counters.PhasedAssemblyMergingTime,
                mContext.Executor, primaryPhaseSets, this::primaryPhasedMerging);
        LOGGER.info("Merged primary phase sets in {}", Counters.PhasedAssemblyMergingTime.formatValue());

        // Secondary phasing
        final List<Set<ExtendedAssembly>> secondaryPhaseSets = Counters.SecondaryPhasingTime.time(
                () -> SecondaryPhasing.run(mergedPhaseSets));
        LOGGER.info("Created {} secondary phase sets in {}", secondaryPhaseSets.size(), Counters.SecondaryPhasingTime.formatValue());

        // Secondary merging
        final List<GappedAssembly> mergedSecondaries = new ArrayList<>();
        Counters.MergeSecondaryTime.time(() ->
        {
            // FIXME: Gapped assemblies
            //for(int i = 0; i < secondaryPhaseSets.size(); i++)
            //    mergedSecondaries.add(createGapped(secondaryPhaseSets.get(i), i));
            for(int i = 0; i < secondaryPhaseSets.size(); i++)
            {
                final Set<ExtendedAssembly> phaseSet = secondaryPhaseSets.get(i);
                int j = 0;
                for(final ExtendedAssembly assembly : phaseSet)
                {
                    final GappedAssembly newAssembly = new GappedAssembly(String.format("Assembly%s-%s", i, j++), List.of(assembly));
                    newAssembly.addErrata(assembly.getAllErrata());
                    for(final Map.Entry<Record, Integer> entry : assembly.getSupport())
                        newAssembly.addEvidenceAt(entry.getKey(), entry.getValue());

                    mergedSecondaries.add(newAssembly);
                }
            }
        });
        LOGGER.info("Merged secondaries in {}", Counters.MergeSecondaryTime.formatValue());

        // Alignment
        final List<AlignedAssembly> aligned = ParallelMapper.mapWithProgress("Alignment", Counters.AlignmentTime,
                mContext.Executor, mergedSecondaries, mContext.Aligner::align);
        LOGGER.info("Created {} alignments in {}", aligned.size(), Counters.AlignmentTime.formatValue());

        // Left sliding (we will find mid-points after calling + de-duping)
        final List<AlignedAssembly> homologised = ParallelMapper.mapWithProgress("Homology Sliding", Counters.HomologyTime,
                mContext.Executor, aligned, mHomologySlider::slideHomology);
        LOGGER.info("Processed homology in {}", Counters.HomologyTime.formatValue());

        // Support scanning
        final SupportScanner supportScanner = new SupportScanner(mContext, Counters.ExtraScannedSupport);
        ParallelMapper.mapWithProgress("Support Scan", Counters.SupportScanTime,
                mContext.Executor, homologised, supportScanner::tryRescanSupport);
        LOGGER.info("Rescanned support, adding {} new reads in {}", Counters.ExtraScannedSupport.formatValue(),
                Counters.SupportScanTime.formatValue());

        // Calling
        final List<VariantCall> variants = Counters.VariantCallingTime.time(() -> new VariantCaller(mContext.Config, mContext.Executor)
                .callVariants(homologised));
        LOGGER.info("Called {} variants in {}", variants.size(), Counters.VariantCallingTime.formatValue());

        // Variant Deduplication
        final VariantDeduplication deduplicator = new VariantDeduplication(mContext, Counters.VariantDeduplicationCounters);
        final List<VariantCall> deduplicated = deduplicator.deduplicate(variants);

        LOGGER.info("{} variants remaining after deduplication", deduplicated.size());
        deduplicated.removeIf(variant -> variant.supportingFragments().isEmpty());
        LOGGER.info("{} variants remaining after removing unsubstantiated", deduplicated.size());

        final long lowQualityVariants = deduplicated.stream()
                .filter(variant -> variant.quality() < mContext.Config.vcfLowQualityThreshold())
                .count();
        LOGGER.info("{} low-quality variants found", lowQualityVariants);

        final long lowSupportVariants = deduplicated.stream()
                .filter(variant -> variant.quality() >= mContext.Config.vcfLowQualityThreshold())
                .filter(variant -> variant.supportingFragments().size() < mContext.Config.minReadsToSupportAssembly())
                .count();
        LOGGER.info("{} low-support variants found (excl low-quality)", lowSupportVariants);

        if(mContext.Config.dropGermline())
        {
            deduplicated.removeIf(VariantCall::isGermline);
            LOGGER.info("{} variants remaining after dropping those with germline support", deduplicated.size());
        }

        if(!mContext.Problems.isEmpty())
        {
            LOGGER.warn("Encountered {} problems", mContext.Problems.size());
            if(mContext.Problems.size() < 50)
            {
                for(final Problem problem : mContext.Problems)
                    LOGGER.warn("{}", problem);
            }
        }
        else
            LOGGER.info("No problems encountered");
        final long endTimeNanos = System.nanoTime();
        LOGGER.info("Completed processing in {}", StringUtils.formatNanos(endTimeNanos - startTimeNanos));

        if(mContext.Config.createHTMLSummaries())
            writeHTMLSummaries(deduplicated);

        writeVCF(deduplicated);

        if(mContext.Config.outputCSVFile() != null)
            writeCSV(deduplicated);

        if(mContext.Config.outputBAMFile() != null)
            writeBAM(deduplicated);

        return deduplicated;
    }

    private void writeHTMLSummaries(final List<VariantCall> variants)
    {
        int summariesWritten = 0;
        for(final VariantCall call : variants)
        {
            if(summariesWritten++ > mContext.Config.maxHTMLSummaries())
            {
                LOGGER.warn("Not writing further HTML summaries -- limit reached. Increase -max_html_summaries to see more.");
                break;
            }

            try
            {
                VariantCallPageGenerator.generatePage(mContext.Config.htmlSummariesFolder(), mContext.ReferenceGenome, mContext.SupportChecker, call);
            }
            catch(final Exception ex)
            {
                LOGGER.error("Failed to generate HTML for {}", call, ex);
            }
        }

        try
        {
            SummaryPageGenerator.generatePage(mContext.Config.htmlSummariesFolder(), Counters, variants);
        }
        catch(final Exception ex)
        {
            LOGGER.error("Failure while generating summary HTML", ex);
        }
    }

    private void writeCSV(final List<VariantCall> variants)
    {
        variants.sort(Comparator.<VariantCall, String>comparing(v -> Objects.requireNonNullElse(v.LeftChromosome, v.RightChromosome))
                .thenComparingInt(v -> v.LeftPosition == 0 ? v.RightPosition : v.LeftPosition));

        CSVWriter.writeCSV(mContext.Config.outputCSVFile(), VariantLine.class, variants.stream()
                .map(variant -> new VariantLine(variant.LeftChromosome, variant.LeftPosition,
                        variant.RightChromosome, variant.RightPosition,
                        variant.LeftDescriptor, variant.RightDescriptor,
                        variant.LeftMappingQuality,
                        variant.germlineSupport(), variant.somaticSupport(),
                        variant.associatedAssemblies().stream().map(asm -> asm.Assembly).collect(Collectors.toList()),
                        variant.Classification,
                        List.of(), csvFilters(variant))));
    }

    private String csvFilters(final VariantCall variant)
    {
        final boolean isLowOverhang = variant.overhang() < mContext.Config.vcfLowOverhangThreshold();
        final boolean isLowQuality = variant.quality() < mContext.Config.vcfLowQualityThreshold();
        final boolean isLowSupport = variant.supportingFragments().size() < mContext.Config.minReadsToSupportAssembly();
        final boolean isLikelyFalse = isLowSupport || (isLowOverhang && variant.discordantSupport() == 0) || isLowQuality;
        final List<String> filters = new ArrayList<>();
        if(variant.isGermline())
            filters.add("GERMLINE");
        if(variant.associatedAssemblies().size() > 1)
            filters.add("MULTIPLE_ASSEMBLIES");
        if(isLowOverhang)
            filters.add("LOW_OVERHANG");
        if(isLowQuality)
            filters.add("LOW_QUALITY");
        if(isLowSupport)
            filters.add("LOW_SUPPORT");
        if(isLikelyFalse)
            filters.add("LIKELY_FALSE");
        return String.join(";", filters);
    }

    private void writeVCF(final List<VariantCall> variants)
    {
        final List<String> sampleNames = variants.stream()
                .flatMap(call -> call.sampleSupport().stream().map(VariantCall.SampleSupport::sampleName))
                .distinct()
                .sorted()
                .collect(Collectors.toList());
        final VCFWriter writer = new VCFWriter(mContext, sampleNames);
        for(final VariantCall call : variants)
        {
            try
            {
                writer.append(call);
            }
            catch(final Exception ex)
            {
                LOGGER.error("Failure while appending to call VCF: {}", call, ex);
            }
        }
        writer.close();
    }

    private void writeBAM(final List<VariantCall> variants)
    {
        final Map<String, List<VariantBAMRecord>> recordsByAssemblyName = new HashMap<>();
        for(final VariantCall call : variants)
        {
            final Set<String> supportFragments = call.sampleSupport().stream()
                    .flatMap(s -> Stream.concat(s.splitReads().stream(), s.discordantReads().stream()))
                    .map(Record::getName)
                    .collect(Collectors.toSet());
            for(final VariantCall.VariantAssembly variantAssembly : call.variantAssemblies())
            {
                final AlignedAssembly assembly = variantAssembly.Assembly;
                final List<VariantBAMRecord> assemblyRecords = recordsByAssemblyName.computeIfAbsent(assembly.Name, __ ->
                {
                    final List<VariantBAMRecord> records = new ArrayList<>();
                    final Set<String> junctions = assembly.getAllErrata(JunctionMetrics.class).stream()
                            .map(junction -> junction.JunctionChromosome + ":" + junction.JunctionPosition
                                    + junction.JunctionDirection.toShortString())
                            .collect(Collectors.toCollection(() -> new TreeSet<>(NaturalSortComparator.INSTANCE)));
                    final Set<String> fragments = new HashSet<>();

                    final List<VariantAssemblyAlignment> alignments = constructAlignments(assembly.getAlignmentBlocks());
                    for(final VariantAssemblyAlignment alignment : alignments)
                    {
                        records.add(new VariantBAMRecord(assembly.Name,
                                assembly.getBasesString().replace('X', 'N'), assembly.getBaseQuality(),
                                alignment.Offset, alignment.Chromosome, alignment.Position, alignment.Inverted,
                                alignment.MapQ,
                                alignment.Cigar,
                                new HashSet<>(), junctions, fragments,
                                alignments));
                    }
                    return records;
                });
                final Set<String> support = new HashSet<>();
                support.addAll(supportFragments);
                support.retainAll(assembly.getSupportFragments());
                assemblyRecords.forEach(r -> r.SourceFragments.addAll(support));
                assemblyRecords.forEach(r -> r.Variants.add(call.compactName()));
            }
        }

        final List<VariantBAMRecord> records = recordsByAssemblyName.values().stream()
                .flatMap(Collection::stream)
                .sorted(NaturalSortComparator.<VariantBAMRecord>of(r -> r.Chromosome).thenComparingInt(r -> r.Position))
                .collect(Collectors.toList());

        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        mContext.ReferenceGenome.refGenomeFile().getSequenceDictionary().getSequences().stream()
                .sorted(NaturalSortComparator.of(SAMSequenceRecord::getSequenceName))
                .map(seq -> new SAMSequenceRecord(seq.getSequenceName(), seq.getSequenceLength()))
                .forEach(header::addSequence);
        final String outputFilename = mContext.Config.outputBAMFile();
        assert outputFilename != null;
        try
        {
            final var writer = new BAMStreamWriter(new FileOutputStream(outputFilename),
                    new FileOutputStream(outputFilename + ".bai"),
                    null, 0,
                    header);

            writer.writeHeader(header);
            for (final VariantBAMRecord record : records)
            {
                final var samRecord = new SAMRecord(header);
                samRecord.setReadName(record.Name);
                samRecord.setReadBases(record.Assembly.getBytes());
                samRecord.setBaseQualities(record.Quality);
                samRecord.setReferenceName(record.Chromosome);
                samRecord.setAlignmentStart(record.Position);
                samRecord.setReadNegativeStrandFlag(record.Inverted);
                samRecord.setMappingQuality(record.MapQ);
                samRecord.setCigarString(record.Cigar);
                samRecord.setAttribute("XV", String.join(",", record.Variants));
                samRecord.setAttribute("XJ", String.join(",", record.Junctions));
                samRecord.setAttribute("XF", String.join(",", record.SourceFragments));
                @Nullable
                final VariantAssemblyAlignment mate = record.getMate();
                if (mate != null)
                {
                    samRecord.setReadPairedFlag(true);
                    samRecord.setFirstOfPairFlag(record.Offset < mate.Offset);
                    samRecord.setSecondOfPairFlag(record.Offset > mate.Offset);
                    samRecord.setMateReferenceName(mate.Chromosome);
                    samRecord.setMateAlignmentStart(mate.Position);
                    samRecord.setMateNegativeStrandFlag(mate.Inverted);
                    samRecord.setAttribute("MC", mate.Cigar);
                }
                if (!record.Supplementaries.isEmpty())
                {
                    final String supplementals = record.Supplementaries.stream()
                            .map(sup -> String.join(",", sup.Chromosome, String.valueOf(sup.Position),
                                    sup.Inverted ? "-" : "+", sup.Cigar, String.valueOf(sup.MapQ), "0") + ";")
                            .collect(Collectors.joining());
                    samRecord.setAttribute("SA", supplementals);
                }

                writer.writeAlignment(samRecord);
            }

            writer.finish(true);
        }
        catch(final Exception ex)
        {
            LOGGER.warn("Failure while writing output BAM", ex);
        }
    }

    private List<VariantAssemblyAlignment> constructAlignments(final List<Alignment> alignments)
    {
        final Set<Alignment> used = Collections.newSetFromMap(new IdentityHashMap<>());

        final List<VariantAssemblyAlignment> variantAlignments = new ArrayList<>();
        for(int i = 0; i < alignments.size(); i++)
        {
            final Alignment alignment = alignments.get(i);
            if(alignment.isUnmapped() || used.contains(alignment))
                continue;
            variantAlignments.add(constructAlignment(alignments, i, used));
        }

        return variantAlignments;
    }

    private VariantAssemblyAlignment constructAlignment(final List<Alignment> alignments, final int startIndex, final Set<Alignment> used)
    {
        final List<CigarElement> elements = new ArrayList<>();

        final Alignment start = alignments.get(startIndex);
        if(start.SequenceStartPosition != 1)
            elements.add(new CigarElement(start.SequenceStartPosition - 1, CigarOperator.S));
        int index = startIndex;
        for(; index < alignments.size(); index++)
        {
            final Alignment current = alignments.get(index);
            used.add(current);
            elements.add(new CigarElement(current.Length, CigarOperator.M));

            if(index + 1 < alignments.size())
            {
                final Alignment next = alignments.get(index + 1);
                if(next.isMapped())
                {
                    if(next.Inverted != current.Inverted || !next.Chromosome.equals(current.Chromosome))
                        break;

                    // Might be a delete
                    final int impliedDelete = !current.Inverted
                            ? next.ReferenceStartPosition - (current.ReferenceStartPosition + current.Length)
                            : current.ReferenceStartPosition - (next.ReferenceStartPosition + next.Length);
                    if(impliedDelete > 0 && impliedDelete <= 100)
                    {
                        elements.add(new CigarElement(impliedDelete, CigarOperator.D));
                        continue;
                    }
                }
            }
            if(index + 2 < alignments.size())
            {
                // Does a future alignment come back to within 100 bases of us?
                int insertSize = alignments.get(index + 1).Length;
                boolean addedInsert = false;
                for(int checkIndex = index + 2; checkIndex < alignments.size(); checkIndex++)
                {
                    final Alignment next = alignments.get(checkIndex);
                    if(next.isUnmapped() || next.Inverted != current.Inverted || !next.Chromosome.equals(current.Chromosome))
                        continue;

                    final int impliedDelete = !current.Inverted
                            ? next.ReferenceStartPosition - (current.ReferenceStartPosition + current.Length)
                            : current.ReferenceStartPosition - (next.ReferenceStartPosition + next.Length);
                    if(impliedDelete >= 0 && impliedDelete < 20)
                    {
                        elements.add(new CigarElement(insertSize, CigarOperator.I));
                        if(impliedDelete > 0)
                            elements.add(new CigarElement(impliedDelete, CigarOperator.D));
                        index = checkIndex - 1;
                        addedInsert = true;
                        break;
                    }

                    insertSize += next.Length;
                }
                if (addedInsert)
                    continue;
            }
            break;
        }

        int softClipSize = 0;
        for(index++; index < alignments.size(); index++)
            softClipSize += alignments.get(index).Length;
        if(softClipSize != 0)
            elements.add(new CigarElement(softClipSize, CigarOperator.S));
        if (start.Inverted)
            Collections.reverse(elements);

        return new VariantAssemblyAlignment(startIndex, start.Chromosome, start.ReferenceStartPosition,
                start.Inverted, start.Quality, new Cigar(elements).toString());
    }

    private static class VariantAssemblyAlignment
    {
        public final int Offset;
        public final String Chromosome;
        public final int Position;
        public final boolean Inverted;
        public final int MapQ;
        public final String Cigar;

        private VariantAssemblyAlignment(final int offset, final String chromosome, final int position,
                final boolean inverted, final int mapQ, final String cigar)
        {
            Offset = offset;
            Chromosome = chromosome;
            Position = position;
            Inverted = inverted;
            MapQ = mapQ;
            Cigar = cigar;
        }
    }

    private static class VariantBAMRecord
    {
        public final String Name;
        public final String Assembly;
        public final byte[] Quality;
        public final int Offset;
        public final String Chromosome;
        public final int Position;
        public final boolean Inverted;
        public final int MapQ;
        public final String Cigar;
        public final Set<String> Variants;
        public final Set<String> Junctions;
        public final Set<String> SourceFragments;
        public final List<VariantAssemblyAlignment> Supplementaries;

        private VariantBAMRecord(final String name, final String assembly, final byte[] quality, final int offset, final String chromosome,
                final int position, final boolean inverted, final int mapQ, final String cigar, final Set<String> variants,
                final Set<String> junctions, final Set<String> sourceFragments, final List<VariantAssemblyAlignment> supplementaries)
        {
            Name = name;
            if (inverted)
            {
                Assembly = SequenceUtil.reverseComplement(assembly);
                final byte[] reversedQuals = Arrays.copyOf(quality, quality.length);
                SequenceUtil.reverseQualities(quality);
                Quality = reversedQuals;
            }
            else
            {
                Assembly = assembly;
                Quality = quality;
            }
            Offset = offset;
            Chromosome = chromosome;
            Position = position;
            Inverted = inverted;
            MapQ = mapQ;
            Cigar = cigar;
            Variants = variants;
            Junctions = junctions;
            SourceFragments = sourceFragments;
            Supplementaries = supplementaries.stream()
                    .filter(s -> s.Offset != Offset)
                    .collect(Collectors.toList());
        }

        @Nullable
        public VariantAssemblyAlignment getMate()
        {
            if(Supplementaries.size() != 1)
                return null;

            return Supplementaries.get(0);
        }
    }

    private Set<ExtendedAssembly> primaryPhasedMerging(final Set<ExtendedAssembly> primaryPhaseSet)
    {
        final Timeout timeout = new Timeout(mContext.Config, TimeUnit.SECONDS.toNanos(5));
        try
        {
            final Set<ExtendedAssembly> result = new HashSet<>(primaryPhaseSet);
            final Set<Pair<ExtendedAssembly, ExtendedAssembly>> checked = new HashSet<>();
            boolean merged = true;
            while(merged)
            {
                merged = false;

                loopHead:
                for(final ExtendedAssembly left : result)
                    for(final ExtendedAssembly right : result)
                    {
                        timeout.checkTimeout();

                        if(left == right)
                            continue;
                        if(!checked.add(Pair.of(left, right)))
                            continue;

                        final int minOverlap = Math.min(30, Math.min(left.getLength(), right.getLength()));
                        @Nullable
                        Integer index = mContext.SupportChecker.AssemblySupport.supportIndex(left, right, minOverlap);
                        if (index != null)
                            index = mContext.SupportChecker.AssemblySupport.bestSupportIndex(left, right, minOverlap);
                        final ExtendedAssembly mergedAssembly;
                        if(index == null)
                        {
                            final ExtendedAssembly flippedRight = right.flipStrand();
                            index = mContext.SupportChecker.AssemblySupport.supportIndex(left, flippedRight, minOverlap);
                            if (index != null)
                                index = mContext.SupportChecker.AssemblySupport.bestSupportIndex(left, right, minOverlap);
                            if(index == null)
                                continue;

                            mergedAssembly = merge(left, flippedRight, index);
                        }
                        else
                            mergedAssembly = merge(left, right, index);

                        result.remove(left);
                        result.remove(right);
                        result.add(mergedAssembly);

                        merged = true;
                        break loopHead;
                    }
            }

            return result;
        }
        catch(final Throwable throwable)
        {
            LOGGER.warn("Failure during phased assembly merging with group of size {}", primaryPhaseSet.size(), throwable);
            LOGGER.warn("{}", RegionOfInterest.tryMerge(
                    primaryPhaseSet.stream()
                            .flatMap(assembly -> assembly.getSupport().stream())
                            .map(Map.Entry::getKey)
                            .filter(record -> !record.isUnmapped())
                            .map(record -> new RegionOfInterest(record.getChromosome(), record.getAlignmentStart(), record.getAlignmentEnd()))
                            .collect(Collectors.toList())
            ));
            return null;
        }
    }

    private ExtendedAssembly merge(final ExtendedAssembly left, final ExtendedAssembly right, final int supportIndex)
    {
        left.markDecompositionStale();
        right.markDecompositionStale();
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final ExtendedAssembly merged = new ExtendedAssembly(left.Name, mergedSequence.getBasesString(), left.Source);
        left.Diagrams.forEach(merged::addDiagrams);

        left.getSupportRecords().forEach(support -> merged.tryAddSupport(mContext.SupportChecker, support));
        right.getSupportRecords().forEach(support -> merged.tryAddSupport(mContext.SupportChecker, support));

        merged.addErrata(left.getAllErrata());
        merged.addErrata(right.getAllErrata());

        return merged;
    }

    private AlignedAssembly merge(final AlignedAssembly left, final AlignedAssembly right, final int supportIndex)
    {
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final ExtendedAssembly merged = new ExtendedAssembly(left.Name, mergedSequence.getBasesString(), left.Source);
        left.Source.Sources.get(0).Diagrams.forEach(merged::addDiagrams);

        final GappedAssembly gapped = new GappedAssembly(merged.Name, List.of(merged));
        reAddSupport(gapped, left);
        reAddSupport(gapped, right);

        return mHomologySlider.slideHomology(mContext.Aligner.align(gapped));
    }

    private void reAddSupport(final SupportedAssembly merged, final SupportedAssembly old)
    {
        final int offset = merged.Assembly.indexOf(old.Assembly);
        for(final Map.Entry<Record, Integer> entry : old.getSupport())
        {
            final Record potentialSupport = entry.getKey();
            if(offset != -1)
            {
                final int oldSupportIndex = entry.getValue();
                if(mContext.SupportChecker.AssemblySupport.supportsAt(merged, potentialSupport, oldSupportIndex + offset))
                {
                    merged.addEvidenceAt(potentialSupport, oldSupportIndex + offset);
                    continue;
                }
            }
            merged.tryAddSupport(mContext.SupportChecker, potentialSupport);
        }
    }

    public List<PrimaryAssembly> flatten(final List<PrimaryAssemblyResult> results)
    {
        return results.stream()
                .flatMap(r -> r.Assemblies.stream())
                .collect(Collectors.toList());
    }

    public List<PrimaryAssembly> consolidateNearbyAssemblies(final List<PrimaryAssembly> results)
    {
        final List<PrimaryAssembly> firstPass = consolidateNearbyAssemblies(results, this::tryMerge);
        return consolidateNearbyAssemblies(firstPass, (left, right) ->
        {
            final Set<String> leftOnly = new HashSet<>(left.getSupportFragments());
            leftOnly.removeAll(right.getSupportFragments());

            final Set<String> rightOnly = new HashSet<>(right.getSupportFragments());
            rightOnly.removeAll(left.getSupportFragments());

            if (leftOnly.size() < 2 && rightOnly.size() < 2)
            {
                final boolean returnRight;
                if(leftOnly.size() == rightOnly.size())
                {
                    if(left.getAverageBaseQuality() == right.getAverageBaseQuality())
                        returnRight = left.getLength() < right.getLength();
                    else
                        returnRight = left.getAverageBaseQuality() < right.getAverageBaseQuality();
                }
                else
                    returnRight = leftOnly.size() < rightOnly.size();

                if(returnRight)
                {
                    right.addErrata(left.getAllErrata());
                    return right;
                }
                else
                {
                    left.addErrata(right.getAllErrata());
                    return left;
                }
            }
            if (leftOnly.size() < 2)
            {
                right.addErrata(left.getAllErrata());
                return right;
            }
            else if (rightOnly.size() < 2)
            {
                left.addErrata(right.getAllErrata());
                return left;
            }
            else
                return null;
        });
    }

    public List<PrimaryAssembly> consolidateNearbyAssemblies(final List<PrimaryAssembly> results,
            final BiFunction<PrimaryAssembly, PrimaryAssembly, PrimaryAssembly> merger)
    {
        results.sort(NaturalSortComparator.<PrimaryAssembly>of(r -> r.AnchorChromosome)
                .thenComparing(r -> r.AnchorPosition));

        final List<PrimaryAssembly> assemblies = new ArrayList<>();
        for(int i = 0; i < results.size(); i++)
        {
            @Nullable
            PrimaryAssembly current = results.get(i);
            if(current == null)
                continue;

            final int maxDedupeDistance = mContext.Config.maxDistanceToDedupeAssemblies();
            final int maxToCheck = Math.min(results.size() - i - 1, maxDedupeDistance * 2);
            for(int j = 0; j < maxToCheck; j++)
            {
                final PrimaryAssembly next = results.get(i + j + 1);
                if(next == null)
                    continue;
                if(!current.AnchorChromosome.equals(next.AnchorChromosome))
                    break;

                final int currentStart = current.AnchorPosition - current.AnchorPositionInAssembly;
                final int nextStart = next.AnchorPosition - next.AnchorPositionInAssembly;
                if(!RangeUtils.overlaps(currentStart - maxDedupeDistance, currentStart + current.Assembly.length() + maxDedupeDistance,
                        nextStart - maxDedupeDistance, nextStart + next.Assembly.length() + maxDedupeDistance))
                    break;

                try
                {
                    @Nullable
                    final PrimaryAssembly merged = merger.apply(current, next);
                    if(merged != null)
                    {
                        current = merged;
                        results.set(i + j + 1, null); // Null out next
                    }
                }
                catch(final Throwable throwable)
                {
                    mContext.Problems.add(new Problem(String.format("Problem merging %s and %s", current.Name, next.Name),
                            throwable, current));
                }
            }

            assemblies.add(current);
        }
        return assemblies;
    }

    @Nullable
    public PrimaryAssembly tryMerge(final PrimaryAssembly left, final PrimaryAssembly right)
    {
        @Nullable
        final Integer mergeIndex = mContext.SupportChecker.AssemblySupport.supportIndex(left, right, 100);
        if(mergeIndex == null)
            return null;

        return merge(left, right, mergeIndex);
    }

    private PrimaryAssembly merge(final PrimaryAssembly left, final PrimaryAssembly right, final int supportIndex)
    {
        final Sequence mergedSequence = SequenceMerger.merge(left, right, supportIndex);

        final var merged = new PrimaryAssembly(left.Name, mergedSequence.getBasesString(),
                "?", 0, 0, left);

        final int leftDelta = supportIndex > 0 ? 0 : -supportIndex;
        for(final Map.Entry<Record, Integer> entry : left.getSupport())
            merged.tryAddSupport(mContext.SupportChecker, entry.getKey(), entry.getValue() + leftDelta);

        final int rightDelta = Math.max(supportIndex, 0);
        for(final Map.Entry<Record, Integer> entry : right.getSupport())
            merged.tryAddSupport(mContext.SupportChecker, entry.getKey(), entry.getValue() + rightDelta);

        merged.addErrata(left.getAllErrata());
        merged.addErrata(right.getAllErrata());

        return merged;
    }

    private List<ExtendedAssembly> order(final Collection<ExtendedAssembly> assemblies)
    {
        // FIXME: Correctly order these
        if(assemblies.size() > 1)
            LOGGER.warn("Found more than 1 assembly ({}) while creating gapped ({})", assemblies.size(),
                    assemblies.stream().map(assembly -> assembly.Name).collect(Collectors.toList()));

        final Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> leftWise = new IdentityHashMap<>();
        final Map<ExtendedAssembly, Map<ExtendedAssembly, Long>> rightWise = new IdentityHashMap<>();
        for (final ExtendedAssembly first : assemblies)
            for (final ExtendedAssembly second : assemblies)
            {
                if (first == second)
                    continue;



            }

        return new ArrayList<>(assemblies);
    }

    public GappedAssembly createGapped(final Collection<ExtendedAssembly> assemblies, final int index)
    {
        final GappedAssembly gappedAssembly = new GappedAssembly("Assembly" + index, order(assemblies));

        for(final ExtendedAssembly assembly : assemblies)
            for(final Record support : assembly.getSupportRecords())
                if(!gappedAssembly.tryAddSupport(mContext.SupportChecker, support))
                    LOGGER.info("Failed to add support for assembly {}: {}", gappedAssembly.Name, support.getName());

        return gappedAssembly;
    }
}
