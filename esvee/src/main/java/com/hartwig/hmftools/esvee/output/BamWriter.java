package com.hartwig.hmftools.esvee.output;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.old.AlignedAssembly;
import com.hartwig.hmftools.esvee.old.Alignment;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.old.VariantAssembly;
import com.hartwig.hmftools.esvee.variant.VariantCall;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.BAMStreamWriter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SequenceUtil;

public class BamWriter
{
    private final BAMStreamWriter mWriter;
    private SAMFileHeader mHeader;

    public BamWriter(final SvConfig config)
    {
        mHeader = null;
        mWriter = initialiseBam(config);
    }

    private BAMStreamWriter initialiseBam(final SvConfig config)
    {
        if(!config.WriteTypes.contains(WriteType.ASSEMBLY_BAM))
            return null;

        mHeader = new SAMFileHeader();
        mHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        RefGenomeSource refGenomeFile = (RefGenomeSource)config.RefGenome;

        refGenomeFile.refGenomeFile().getSequenceDictionary().getSequences().stream()
                .sorted(NaturalSortComparator.of(SAMSequenceRecord::getSequenceName))
                .map(seq -> new SAMSequenceRecord(seq.getSequenceName(), seq.getSequenceLength()))
                .forEach(mHeader::addSequence);

        String outputFilename = config.outputFilename(WriteType.ASSEMBLY_BAM);

        try
        {
            final BAMStreamWriter writer = new BAMStreamWriter(new FileOutputStream(outputFilename),
                    new FileOutputStream(outputFilename + ".bai"),
                    null, 0, mHeader);

            writer.writeHeader(mHeader);

            return writer;
        }
        catch(final Exception e)
        {
            SV_LOGGER.warn("failure while writing output BAM", e.toString());
            return null;
        }
    }

    public void close()
    {
        if(mWriter != null)
            mWriter.finish(true);
    }

    public void writeVariantAssemblyBamRecords(final List<VariantCall> variants)
    {
        if(mWriter == null)
            return;

        Map<String, List<VariantBAMRecord>> recordsByAssemblyName = new HashMap<>();
        
        for(VariantCall call : variants)
        {
            Set<String> supportFragments = call.sampleSupport().stream()
                    .flatMap(s -> Stream.concat(s.splitReads().stream(), s.discordantReads().stream()))
                    .map(Read::getName)
                    .collect(Collectors.toSet());

            for(VariantAssembly variantAssembly : call.variantAssemblies())
            {
                final AlignedAssembly assembly = variantAssembly.Assembly;
                final List<VariantBAMRecord> assemblyRecords = recordsByAssemblyName.computeIfAbsent(assembly.Name, __ ->
                {
                    final List<VariantBAMRecord> records = new ArrayList<>();

                    /*
                    final Set<String> junctions = assembly.getAllErrata(JunctionMetrics.class).stream()
                            .map(junction -> junction.JunctionChromosome + ":" + junction.JunctionPosition
                                    + junction.JunctionDirection.toShortString())
                            .collect(Collectors.toCollection(() -> new TreeSet<>(NaturalSortComparator.INSTANCE)));
                    */

                    Set<String> junctions = Sets.newHashSet();

                    final Set<String> fragments = new HashSet<>();

                    final List<VariantAssemblyAlignment> alignments = constructAlignments(assembly.getAlignmentBlocks());

                    for(VariantAssemblyAlignment alignment : alignments)
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
                support.retainAll(assembly.getSupportReadNames());
                assemblyRecords.forEach(r -> r.SourceFragments.addAll(support));
                assemblyRecords.forEach(r -> r.Variants.add(call.compactName()));
            }
        }

        final List<VariantBAMRecord> records = recordsByAssemblyName.values().stream()
                .flatMap(Collection::stream)
                .sorted(NaturalSortComparator.<VariantBAMRecord>of(r -> r.Chromosome).thenComparingInt(r -> r.Position))
                .collect(Collectors.toList());

        try
        {
            for(VariantBAMRecord record : records)
            {
                final var samRecord = new SAMRecord(mHeader);
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
                if(mate != null)
                {
                    samRecord.setReadPairedFlag(true);
                    samRecord.setFirstOfPairFlag(record.Offset < mate.Offset);
                    samRecord.setSecondOfPairFlag(record.Offset > mate.Offset);
                    samRecord.setMateReferenceName(mate.Chromosome);
                    samRecord.setMateAlignmentStart(mate.Position);
                    samRecord.setMateNegativeStrandFlag(mate.Inverted);
                    samRecord.setAttribute("MC", mate.Cigar);
                }
                
                if(!record.Supplementaries.isEmpty())
                {
                    final String supplementals = record.Supplementaries.stream()
                            .map(sup -> String.join(",", sup.Chromosome, String.valueOf(sup.Position),
                                    sup.Inverted ? "-" : "+", sup.Cigar, String.valueOf(sup.MapQ), "0") + ";")
                            .collect(Collectors.joining());
                    samRecord.setAttribute("SA", supplementals);
                }

                mWriter.writeAlignment(samRecord);
            }
        }
        catch(final Exception e)
        {
            SV_LOGGER.warn("failure while writing output BAM", e.toString());
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

        private VariantBAMRecord(
                final String name, final String assembly, final byte[] quality, final int offset, final String chromosome,
                final int position, final boolean inverted, final int mapQ, final String cigar, final Set<String> variants,
                final Set<String> junctions, final Set<String> sourceFragments, final List<VariantAssemblyAlignment> supplementaries)
        {
            Name = name;
            if(inverted)
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
            Supplementaries = supplementaries.stream().filter(s -> s.Offset != Offset).collect(Collectors.toList());
        }

        @Nullable
        public VariantAssemblyAlignment getMate()
        {
            if(Supplementaries.size() != 1)
                return null;

            return Supplementaries.get(0);
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
                if(addedInsert)
                    continue;
            }
            break;
        }

        int softClipSize = 0;
        for(index++; index < alignments.size(); index++)
            softClipSize += alignments.get(index).Length;
        if(softClipSize != 0)
            elements.add(new CigarElement(softClipSize, CigarOperator.S));
        if(start.Inverted)
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
}
